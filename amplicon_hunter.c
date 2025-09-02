#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <libgen.h>
#include <math.h>
#include <immintrin.h>  // Added for AVX2
#ifdef __linux__
#include <malloc.h>     // For malloc_trim on Linux
#endif

#define VERSION "2.0.0"  // IUPAC-aware + chunked reading version
#define MAX_LINE_LENGTH 65536
#define MAX_PATH_LENGTH 4096
#define MAX_SEQ_LENGTH 100000000
#define DEFAULT_THREADS 8
#define MAX_PRIMER_LENGTH 100
#define BATCH_SIZE 1000
#define DEFAULT_BATCH_SIZE_MB 500
#define AMPLICON_BUFFER_SIZE 100  // Reduced to write more frequently and use less memory
#define MMAP_CHUNK_SIZE (16 * 1024 * 1024)  // Process files in 16MB chunks for better memory control

// DNA base bit representation for IUPAC matching
#define A_BIT 0x01
#define C_BIT 0x02
#define G_BIT 0x04
#define T_BIT 0x08

// 2-bit encoding for DNA storage
#define DNA_A 0
#define DNA_C 1
#define DNA_G 2
#define DNA_T 3

// Memory-mapped file structure
typedef struct {
    void *data;
    size_t size;
    int fd;
} MappedFile;

// Streaming 2bit file structure - CHUNKED READING
typedef struct {
    FILE *fp;
    uint32_t num_sequences;
    uint32_t sequences_read;
    uint8_t *buffer;
    size_t buffer_size;
    size_t buffer_pos;
    size_t buffer_valid;
    uint8_t *partial;  // For incomplete records at chunk boundaries
    size_t partial_size;
} StreamingTwoBitFile;

// Single sequence from 2bit file
typedef struct {
    char *header;
    uint32_t length;
    uint32_t packed_length;
    uint8_t *packed_data;  // Now owned by this struct, must be freed
} TwoBitSequence;

// Quality file structure (binary format) - REMOVED FASTQ SUPPORT
typedef struct {
    uint32_t num_sequences;
    uint32_t *lengths;
    uint8_t **data;
    MappedFile mapped;
} QualityFile;

// Batch file structure for compression
typedef struct {
    FILE *seq_fp;
    char seq_path[MAX_PATH_LENGTH];
    size_t current_size;
    uint32_t seq_count;
    int batch_number;
} BatchFile;

// Primer match structure
typedef struct {
    int type;  // 0=F, 1=F_RC, 2=R, 3=R_RC
    uint32_t start;
    uint32_t end;
    char primer_seq[MAX_PRIMER_LENGTH];
} PrimerMatch;

// Amplicon structure
typedef struct {
    char *sequence;
    char *source_file;
    char *seq_id;
    uint32_t start;
    uint32_t end;
    char orientation[3];
    char *forward_barcode;
    char *reverse_barcode;
    char forward_primer_seq[MAX_PRIMER_LENGTH];
    char reverse_primer_seq[MAX_PRIMER_LENGTH];
} Amplicon;

// Temporary file info for amplicon writing
typedef struct {
    char filepath[MAX_PATH_LENGTH];
    FILE *fp;
    int count;
} TempAmpliconFile;

// Thread data for compression with batching
typedef struct {
    int thread_id;
    char **files;
    int num_files;
    int start_idx;
    int end_idx;
    char *output_dir;
    FILE *list_fp;
    pthread_mutex_t *list_mutex;
    size_t batch_size_bytes;
} CompressThreadData;

// Thread data for amplicon finding with periodic writing
typedef struct {
    int thread_id;
    char **seq_files;
    int num_files;
    int start_idx;
    int end_idx;
    
    // Primer info - including precomputed RC masks
    const uint8_t *forward_mask;
    const uint8_t *reverse_mask;
    const uint8_t *forward_rc_mask;
    const uint8_t *reverse_rc_mask;
    int forward_len;
    int reverse_len;
    
    // Parameters
    int mismatches;
    int clamp;
    int min_length;
    int max_length;
    int fb_len;
    int rb_len;
    bool trim_primers;
    bool include_offtarget;
    
    // Results with periodic writing
    Amplicon *amplicons;
    int num_amplicons;
    int capacity;
    
    // Temporary file for writing amplicons
    TempAmpliconFile temp_file;
    char *temp_dir;
    
    // Statistics
    int total_amplicons_written;
} AmpliconThreadData;

// Global lookup tables
static uint8_t CHAR_TO_MASK[256];
static uint8_t CHAR_TO_2BIT[256];
static char BIT2_TO_BASE[4] = {'A', 'C', 'G', 'T'};

// Initialize lookup tables with IUPAC-aware 2-bit encoding
void init_lookup_tables() {
    memset(CHAR_TO_MASK, 0, 256);
    memset(CHAR_TO_2BIT, 0, 256);
    
    // Basic nucleotides for bit masks
    CHAR_TO_MASK['A'] = CHAR_TO_MASK['a'] = A_BIT;
    CHAR_TO_MASK['C'] = CHAR_TO_MASK['c'] = C_BIT;
    CHAR_TO_MASK['G'] = CHAR_TO_MASK['g'] = G_BIT;
    CHAR_TO_MASK['T'] = CHAR_TO_MASK['t'] = T_BIT;
    CHAR_TO_MASK['U'] = CHAR_TO_MASK['u'] = T_BIT;
    
    // IUPAC ambiguity codes for matching
    CHAR_TO_MASK['R'] = CHAR_TO_MASK['r'] = A_BIT | G_BIT;
    CHAR_TO_MASK['Y'] = CHAR_TO_MASK['y'] = C_BIT | T_BIT;
    CHAR_TO_MASK['S'] = CHAR_TO_MASK['s'] = G_BIT | C_BIT;
    CHAR_TO_MASK['W'] = CHAR_TO_MASK['w'] = A_BIT | T_BIT;
    CHAR_TO_MASK['K'] = CHAR_TO_MASK['k'] = G_BIT | T_BIT;
    CHAR_TO_MASK['M'] = CHAR_TO_MASK['m'] = A_BIT | C_BIT;
    CHAR_TO_MASK['B'] = CHAR_TO_MASK['b'] = C_BIT | G_BIT | T_BIT;
    CHAR_TO_MASK['D'] = CHAR_TO_MASK['d'] = A_BIT | G_BIT | T_BIT;
    CHAR_TO_MASK['H'] = CHAR_TO_MASK['h'] = A_BIT | C_BIT | T_BIT;
    CHAR_TO_MASK['V'] = CHAR_TO_MASK['v'] = A_BIT | C_BIT | G_BIT;
    CHAR_TO_MASK['N'] = CHAR_TO_MASK['n'] = A_BIT | C_BIT | G_BIT | T_BIT;
    
    // 2-bit encoding - default to A for unknown
    for (int i = 0; i < 256; i++) {
        CHAR_TO_2BIT[i] = DNA_A;  // Default to A
    }
    
    // Standard bases
    CHAR_TO_2BIT['A'] = CHAR_TO_2BIT['a'] = DNA_A;
    CHAR_TO_2BIT['C'] = CHAR_TO_2BIT['c'] = DNA_C;
    CHAR_TO_2BIT['G'] = CHAR_TO_2BIT['g'] = DNA_G;
    CHAR_TO_2BIT['T'] = CHAR_TO_2BIT['t'] = DNA_T;
    CHAR_TO_2BIT['U'] = CHAR_TO_2BIT['u'] = DNA_T;
    
    // IUPAC codes - use lexicographically first valid base
    CHAR_TO_2BIT['W'] = CHAR_TO_2BIT['w'] = DNA_A;  // A or T -> A
    CHAR_TO_2BIT['S'] = CHAR_TO_2BIT['s'] = DNA_C;  // C or G -> C
    CHAR_TO_2BIT['M'] = CHAR_TO_2BIT['m'] = DNA_A;  // A or C -> A
    CHAR_TO_2BIT['K'] = CHAR_TO_2BIT['k'] = DNA_G;  // G or T -> G
    CHAR_TO_2BIT['R'] = CHAR_TO_2BIT['r'] = DNA_A;  // A or G -> A
    CHAR_TO_2BIT['Y'] = CHAR_TO_2BIT['y'] = DNA_C;  // C or T -> C
    CHAR_TO_2BIT['B'] = CHAR_TO_2BIT['b'] = DNA_C;  // C,G,T -> C
    CHAR_TO_2BIT['D'] = CHAR_TO_2BIT['d'] = DNA_A;  // A,G,T -> A
    CHAR_TO_2BIT['H'] = CHAR_TO_2BIT['h'] = DNA_A;  // A,C,T -> A
    CHAR_TO_2BIT['V'] = CHAR_TO_2BIT['v'] = DNA_A;  // A,C,G -> A
    CHAR_TO_2BIT['N'] = CHAR_TO_2BIT['n'] = DNA_A;  // All four -> A
}

// Convert primer string to bit mask array
uint8_t* primer_to_mask(const char * restrict primer, int * restrict len) {
    *len = strlen(primer);
    uint8_t *mask = (uint8_t*)malloc(*len);
    
    for (int i = 0; i < *len; i++) {
        mask[i] = CHAR_TO_MASK[(unsigned char)primer[i]];
        if (mask[i] == 0) {
            mask[i] = A_BIT | C_BIT | G_BIT | T_BIT;  // Treat unknown as N
        }
    }
    
    return mask;
}

// Reverse complement a bit mask
uint8_t* reverse_complement_mask(const uint8_t * restrict mask, int len) {
    uint8_t *rc = (uint8_t*)malloc(len);
    
    for (int i = 0; i < len; i++) {
        uint8_t orig = mask[len - 1 - i];
        uint8_t comp = 0;
        
        if (orig & A_BIT) comp |= T_BIT;
        if (orig & T_BIT) comp |= A_BIT;
        if (orig & G_BIT) comp |= C_BIT;
        if (orig & C_BIT) comp |= G_BIT;
        
        rc[i] = comp;
    }
    
    return rc;
}

// Reverse complement a DNA string
char* reverse_complement_string(const char * restrict seq, int len) {
    char *rc = (char*)malloc(len + 1);
    
    for (int i = 0; i < len; i++) {
        char c = toupper(seq[len - 1 - i]);
        switch(c) {
            case 'A': rc[i] = 'T'; break;
            case 'T': rc[i] = 'A'; break;
            case 'U': rc[i] = 'A'; break;
            case 'G': rc[i] = 'C'; break;
            case 'C': rc[i] = 'G'; break;
            default: rc[i] = 'N'; break;
        }
    }
    rc[len] = '\0';
    
    return rc;
}

// Process spaces in sequence ID
void process_spaces_in_id(char * restrict id) {
    for (int i = 0; id[i]; i++) {
        if (id[i] == ' ' || id[i] == '\t') id[i] = '_';
    }
}

// Check if file has FASTA extension (NO FASTQ)
bool is_sequence_file(const char * restrict filename) {
    const char *ext = strrchr(filename, '.');
    if (!ext) return false;
    
    // Convert to lowercase for comparison
    char lower_ext[10];
    int i;
    for (i = 0; ext[i] && i < 9; i++) {
        lower_ext[i] = tolower(ext[i]);
    }
    lower_ext[i] = '\0';
    
    return (strcmp(lower_ext, ".fa") == 0 ||
            strcmp(lower_ext, ".fasta") == 0 ||
            strcmp(lower_ext, ".fna") == 0 ||
            strcmp(lower_ext, ".ffn") == 0 ||
            strcmp(lower_ext, ".faa") == 0 ||
            strcmp(lower_ext, ".frn") == 0);
}

// Recursively find all sequence files in directory
void find_sequence_files(const char * restrict dir_path, char *** restrict files, 
                        int * restrict num_files, int * restrict capacity) {
    DIR *dir = opendir(dir_path);
    if (!dir) return;
    
    struct dirent *entry;
    struct stat st;
    char full_path[MAX_PATH_LENGTH];
    
    while ((entry = readdir(dir)) != NULL) {
        if (entry->d_name[0] == '.') continue;
        
        snprintf(full_path, sizeof(full_path), "%s/%s", dir_path, entry->d_name);
        
#ifdef _DIRENT_HAVE_D_TYPE
        if (entry->d_type == DT_DIR) {
            find_sequence_files(full_path, files, num_files, capacity);
        } else if (entry->d_type == DT_REG || entry->d_type == DT_UNKNOWN) {
            if (entry->d_type == DT_UNKNOWN) {
                if (stat(full_path, &st) != 0) continue;
                if (!S_ISREG(st.st_mode)) continue;
            }
            if (is_sequence_file(entry->d_name)) {
                if (*num_files >= *capacity) {
                    *capacity *= 2;
                    *files = (char**)realloc(*files, *capacity * sizeof(char*));
                }
                (*files)[(*num_files)++] = strdup(full_path);
            }
        }
#else
        if (stat(full_path, &st) == 0) {
            if (S_ISDIR(st.st_mode)) {
                find_sequence_files(full_path, files, num_files, capacity);
            } else if (S_ISREG(st.st_mode) && is_sequence_file(entry->d_name)) {
                if (*num_files >= *capacity) {
                    *capacity *= 2;
                    *files = (char**)realloc(*files, *capacity * sizeof(char*));
                }
                (*files)[(*num_files)++] = strdup(full_path);
            }
        }
#endif
    }
    
    closedir(dir);
}

// Write a sequence to a batch file (FASTA ONLY) with IUPAC-aware encoding
size_t write_sequence_to_batch(BatchFile * restrict batch, const char * restrict header, 
                               const char * restrict sequence, uint32_t seq_len) {
    size_t bytes_written = 0;
    
    // Process header
    char processed_header[MAX_LINE_LENGTH];
    strcpy(processed_header, header);
    process_spaces_in_id(processed_header);
    
    // Write header
    uint32_t header_len = strlen(processed_header);
    fwrite(&header_len, sizeof(uint32_t), 1, batch->seq_fp);
    fwrite(processed_header, 1, header_len, batch->seq_fp);
    bytes_written += sizeof(uint32_t) + header_len;
    
    // Write sequence length
    fwrite(&seq_len, sizeof(uint32_t), 1, batch->seq_fp);
    bytes_written += sizeof(uint32_t);
    
    // Pack and write sequence using IUPAC-aware encoding
    uint32_t packed_len = (seq_len + 3) / 4;
    uint8_t *packed = (uint8_t*)calloc(packed_len, 1);
    
    for (uint32_t i = 0; i < seq_len; i++) {
        uint32_t byte_idx = i / 4;
        uint32_t bit_pos = 6 - 2 * (i % 4);
        uint8_t value = CHAR_TO_2BIT[(unsigned char)toupper(sequence[i])];
        packed[byte_idx] |= (value << bit_pos);
    }
    
    fwrite(&packed_len, sizeof(uint32_t), 1, batch->seq_fp);
    fwrite(packed, 1, packed_len, batch->seq_fp);
    free(packed);
    bytes_written += sizeof(uint32_t) + packed_len;
    
    batch->seq_count++;
    batch->current_size += bytes_written;
    
    return bytes_written;
}

// Finalize a batch file
void finalize_batch(BatchFile * restrict batch) {
    if (batch->seq_fp) {
        // Update sequence count at beginning of file
        fseek(batch->seq_fp, 0, SEEK_SET);
        fwrite(&batch->seq_count, sizeof(uint32_t), 1, batch->seq_fp);
        fclose(batch->seq_fp);
        batch->seq_fp = NULL;
    }
}

// Start a new batch file
void start_new_batch(BatchFile * restrict batch, const char * restrict output_dir, 
                    int thread_id) {
    if (batch->seq_fp) {
        finalize_batch(batch);
    }
    
    batch->batch_number++;
    snprintf(batch->seq_path, sizeof(batch->seq_path), "%s/batch_t%d_b%d.2bit", 
             output_dir, thread_id, batch->batch_number);
    
    batch->seq_fp = fopen(batch->seq_path, "wb");
    batch->seq_count = 0;
    batch->current_size = 0;
    
    // Write placeholder for sequence count
    fwrite(&batch->seq_count, sizeof(uint32_t), 1, batch->seq_fp);
}

// Compress files with batching (FASTA ONLY) - STREAMING VERSION
void compress_files_batched(const char * restrict input_file, BatchFile * restrict batch, 
                           const char * restrict output_dir, int thread_id, size_t batch_size_bytes,
                           FILE * restrict list_fp, pthread_mutex_t * restrict list_mutex) {
    // Quick file type check
    int fd = open(input_file, O_RDONLY);
    if (fd < 0) return;
    
    char first_char;
    ssize_t n = read(fd, &first_char, 1);
    close(fd);
    if (n <= 0 || first_char != '>') return;  // FASTA only
    
    FILE *in = fopen(input_file, "r");
    if (!in) return;
    
    char line[MAX_LINE_LENGTH];
    char *sequence = NULL;
    size_t seq_cap = 65536;  // Start with 64KB
    sequence = (char*)malloc(seq_cap);
    uint32_t seq_pos = 0;
    
    // FASTA processing
    char *current_header = NULL;
    
    while (fgets(line, MAX_LINE_LENGTH, in)) {
        if (line[0] == '>') {
            if (seq_pos > 0 && current_header) {
                // Check if we need a new batch
                if (batch->seq_fp == NULL || 
                    batch->current_size + seq_pos + strlen(current_header) + 100 > batch_size_bytes) {
                    if (batch->seq_fp) {
                        finalize_batch(batch);
                        pthread_mutex_lock(list_mutex);
                        fprintf(list_fp, "%s\n", batch->seq_path);
                        pthread_mutex_unlock(list_mutex);
                    }
                    start_new_batch(batch, output_dir, thread_id);
                }
                
                write_sequence_to_batch(batch, current_header, sequence, seq_pos);
                seq_pos = 0;
            }
            
            if (current_header) {
                free(current_header);
                current_header = NULL;
            }
            current_header = strdup(line + 1);  // Skip '>'
            current_header[strcspn(current_header, "\n\r")] = '\0';
        } else {
            int len = strcspn(line, "\n\r");
            // Grow buffer if needed
            if (seq_pos + len >= seq_cap) {
                seq_cap = seq_cap * 2;
                if (seq_cap > MAX_SEQ_LENGTH) seq_cap = MAX_SEQ_LENGTH;
                sequence = (char*)realloc(sequence, seq_cap);
            }
            if (seq_pos + len < seq_cap) {
                memcpy(sequence + seq_pos, line, len);
                seq_pos += len;
            }
        }
    }
    
    // Process last sequence
    if (seq_pos > 0 && current_header) {
        if (batch->seq_fp == NULL || 
            batch->current_size + seq_pos + strlen(current_header) + 100 > batch_size_bytes) {
            if (batch->seq_fp) {
                finalize_batch(batch);
                pthread_mutex_lock(list_mutex);
                fprintf(list_fp, "%s\n", batch->seq_path);
                pthread_mutex_unlock(list_mutex);
            }
            start_new_batch(batch, output_dir, thread_id);
        }
        write_sequence_to_batch(batch, current_header, sequence, seq_pos);
    }
    
    if (current_header) free(current_header);
    free(sequence);
    fclose(in);
}

// Worker thread for compression with batching
void* compress_worker(void *arg) {
    CompressThreadData *data = (CompressThreadData*)arg;
    
    BatchFile batch = {0};
    batch.batch_number = 0;
    
    for (int i = data->start_idx; i < data->end_idx; i++) {
        compress_files_batched(data->files[i], &batch, data->output_dir, 
                             data->thread_id, data->batch_size_bytes,
                             data->list_fp, data->list_mutex);
    }
    
    // Finalize last batch
    if (batch.seq_fp) {
        finalize_batch(&batch);
        pthread_mutex_lock(data->list_mutex);
        fprintf(data->list_fp, "%s\n", batch.seq_path);
        pthread_mutex_unlock(data->list_mutex);
    }
    
    return NULL;
}

// Memory map a file (STREAMING MEMORY MANAGEMENT)
MappedFile map_file(const char * restrict filename) {
    MappedFile mapped = {0};
    mapped.fd = open(filename, O_RDONLY);
    if (mapped.fd == -1) return mapped;
    
    struct stat sb;
    if (fstat(mapped.fd, &sb) == -1) {
        close(mapped.fd);
        mapped.fd = -1;
        return mapped;
    }
    
    mapped.size = sb.st_size;
    // Use MAP_PRIVATE and avoid reserving swap space
    int mmap_flags = MAP_PRIVATE;
#ifdef MAP_NORESERVE
    mmap_flags |= MAP_NORESERVE;  // Don't reserve swap space
#endif
    mapped.data = mmap(NULL, mapped.size, PROT_READ, mmap_flags, mapped.fd, 0);
    
    if (mapped.data == MAP_FAILED) {
        close(mapped.fd);
        mapped.fd = -1;
        mapped.data = NULL;
    } else {
#ifdef __linux__
        // Tell kernel to read ahead sequentially
        posix_madvise(mapped.data, mapped.size, POSIX_MADV_SEQUENTIAL);
#endif
    }
    
    return mapped;
}

// Unmap file (AGGRESSIVE MEMORY RELEASE)
void unmap_file(MappedFile * restrict mapped) {
    if (mapped->data && mapped->data != MAP_FAILED) {
#ifdef __linux__
        // Aggressively tell kernel to drop all pages NOW
        posix_madvise(mapped->data, mapped->size, POSIX_MADV_DONTNEED);
#endif
        munmap(mapped->data, mapped->size);
        mapped->data = NULL;
    }
    if (mapped->fd != -1) {
        close(mapped->fd);
        mapped->fd = -1;
    }
    // Force memory release to OS
#ifdef __linux__
    malloc_trim(0);
#endif
}

// Open a streaming 2bit file with chunked reading
StreamingTwoBitFile open_streaming_2bit(const char * restrict filename) {
    StreamingTwoBitFile file = {0};
    
    file.fp = fopen(filename, "rb");
    if (!file.fp) return file;
    
    // Read number of sequences
    if (fread(&file.num_sequences, sizeof(uint32_t), 1, file.fp) != 1) {
        fclose(file.fp);
        file.fp = NULL;
        return file;
    }
    
    file.sequences_read = 0;
    file.buffer_size = MMAP_CHUNK_SIZE;  // 16MB chunks
    file.buffer = (uint8_t*)malloc(file.buffer_size);
    file.buffer_pos = 0;
    file.buffer_valid = 0;
    file.partial = NULL;
    file.partial_size = 0;
    
    return file;
}

// Helper to ensure we have at least n bytes available in buffer
static bool ensure_bytes(StreamingTwoBitFile * restrict file, size_t n) {
    size_t available = file->buffer_valid - file->buffer_pos;
    
    if (available >= n) {
        return true;
    }
    
    // First, handle any partial data from previous read
    if (file->partial && file->partial_size > 0) {
        // Move partial data to beginning of buffer
        memcpy(file->buffer, file->partial, file->partial_size);
        file->buffer_valid = file->partial_size;
        free(file->partial);
        file->partial = NULL;
        file->partial_size = 0;
    } else {
        // Save remaining data as partial for next read
        if (available > 0) {
            file->partial = (uint8_t*)malloc(available);
            memcpy(file->partial, file->buffer + file->buffer_pos, available);
            file->partial_size = available;
        }
        file->buffer_valid = file->partial_size;
        if (file->partial_size > 0) {
            memcpy(file->buffer, file->partial, file->partial_size);
            free(file->partial);
            file->partial = NULL;
            file->partial_size = 0;
        }
    }
    
    // Read more data from file
    size_t to_read = file->buffer_size - file->buffer_valid;
    size_t bytes_read = fread(file->buffer + file->buffer_valid, 1, to_read, file->fp);
    file->buffer_valid += bytes_read;
    file->buffer_pos = 0;
    
    return file->buffer_valid >= n;
}

// Read next sequence from streaming 2bit file with chunked reading
bool read_next_2bit_sequence(StreamingTwoBitFile * restrict file, TwoBitSequence * restrict seq) {
    if (!file->fp || file->sequences_read >= file->num_sequences) {
        return false;
    }
    
    // Read header length
    if (!ensure_bytes(file, sizeof(uint32_t))) {
        return false;
    }
    uint32_t header_len;
    memcpy(&header_len, file->buffer + file->buffer_pos, sizeof(uint32_t));
    file->buffer_pos += sizeof(uint32_t);
    
    // Read header
    if (!ensure_bytes(file, header_len)) {
        return false;
    }
    seq->header = (char*)malloc(header_len + 1);
    memcpy(seq->header, file->buffer + file->buffer_pos, header_len);
    seq->header[header_len] = '\0';
    file->buffer_pos += header_len;
    
    // Read sequence length
    if (!ensure_bytes(file, sizeof(uint32_t))) {
        free(seq->header);
        return false;
    }
    memcpy(&seq->length, file->buffer + file->buffer_pos, sizeof(uint32_t));
    file->buffer_pos += sizeof(uint32_t);
    
    // Read packed length
    if (!ensure_bytes(file, sizeof(uint32_t))) {
        free(seq->header);
        return false;
    }
    memcpy(&seq->packed_length, file->buffer + file->buffer_pos, sizeof(uint32_t));
    file->buffer_pos += sizeof(uint32_t);
    
    // Read packed data
    if (!ensure_bytes(file, seq->packed_length)) {
        free(seq->header);
        return false;
    }
    seq->packed_data = (uint8_t*)malloc(seq->packed_length);
    memcpy(seq->packed_data, file->buffer + file->buffer_pos, seq->packed_length);
    file->buffer_pos += seq->packed_length;
    
    file->sequences_read++;
    return true;
}

// Free a 2bit sequence (now also frees packed_data)
void free_2bit_sequence(TwoBitSequence * restrict seq) {
    if (seq->header) {
        free(seq->header);
        seq->header = NULL;
    }
    if (seq->packed_data) {
        free(seq->packed_data);
        seq->packed_data = NULL;
    }
}

// Close streaming 2bit file
void close_streaming_2bit(StreamingTwoBitFile * restrict file) {
    if (file->fp) {
        fclose(file->fp);
        file->fp = NULL;
    }
    if (file->buffer) {
        free(file->buffer);
        file->buffer = NULL;
    }
    if (file->partial) {
        free(file->partial);
        file->partial = NULL;
    }
    memset(file, 0, sizeof(StreamingTwoBitFile));
}

// Unpack 2bit sequence
void unpack_sequence(const uint8_t * restrict packed, uint32_t seq_length, char * restrict output) {
    for (uint32_t i = 0; i < seq_length; i++) {
        uint32_t byte_idx = i / 4;
        uint32_t bit_pos = 6 - 2 * (i % 4);
        uint8_t value = (packed[byte_idx] >> bit_pos) & 0x3;
        output[i] = BIT2_TO_BASE[value];
    }
    output[seq_length] = '\0';
}

// SIMD OPTIMIZATION FUNCTIONS START HERE

// Build one-byte mask per base using CHAR_TO_MASK table
static inline void build_seq_mask(const char *s, uint8_t *out, int n) {
    for (int i = 0; i < n; i++) {
        out[i] = CHAR_TO_MASK[(unsigned char)s[i]];
    }
}

// Count positions where (seq_mask[i] & primer_mask[i]) == 0 using AVX2
static inline int count_mismatches_avx2(const uint8_t *seq_mask,
                                        const uint8_t *primer_mask, int n) {
    int mis = 0, i = 0;
    const __m256i z = _mm256_setzero_si256();
    
    for (; i + 32 <= n; i += 32) {
        __m256i a = _mm256_loadu_si256((const __m256i *)(seq_mask + i));
        __m256i b = _mm256_loadu_si256((const __m256i *)(primer_mask + i));
        __m256i c = _mm256_and_si256(a, b);
        __m256i is_zero = _mm256_cmpeq_epi8(c, z);
        unsigned mask = (unsigned)_mm256_movemask_epi8(is_zero);
        mis += __builtin_popcount(mask);
    }
    
    for (; i < n; i++) {
        mis += ((seq_mask[i] & primer_mask[i]) == 0);
    }
    
    return mis;
}

// SIMD version with clamp logic handled on edges
static inline bool matches_with_mismatches_simd(const uint8_t *seq_mask,
                                                const uint8_t *primer_mask,
                                                int len, int max_mis, int clamp,
                                                bool check_rc_clamp) {
    // Enforce clamp exactly on the constrained end
    if (clamp > 0) {
        if (check_rc_clamp) {                    // clamp on 5' (RC path)
            for (int i = 0; i < clamp; i++) {
                if ((primer_mask[i] & seq_mask[i]) == 0) return false;
            }
        } else {                                  // clamp on 3'
            for (int i = len - clamp; i < len; i++) {
                if ((primer_mask[i] & seq_mask[i]) == 0) return false;
            }
        }
    }
    
    // SIMD on the middle segment
    int L = check_rc_clamp ? clamp : 0;
    int R = check_rc_clamp ? len   : len - clamp;
    int mis = 0;
    
    if (R > L) {
        mis += count_mismatches_avx2(seq_mask + L, primer_mask + L, R - L);
        if (mis > max_mis) return false;
    }
    
    return true;
}

// SIMD OPTIMIZATION FUNCTIONS END HERE

// Check if sequence matches primer with mismatches (KEPT FOR FALLBACK/COMPATIBILITY)
bool matches_with_mismatches(const char * restrict seq, const uint8_t * restrict primer_mask, 
                             int len, int max_mismatches, int clamp, bool check_rc_clamp) {
    int mismatches = 0;
    
    for (int i = 0; i < len; i++) {
        uint8_t seq_mask = CHAR_TO_MASK[(unsigned char)seq[i]];
        
        if ((primer_mask[i] & seq_mask) == 0) {
            mismatches++;
            
            // Check if mismatch is in clamp region
            if (check_rc_clamp) {
                if (i < clamp) return false;  // RC: clamp is at start (5' end after RC)
            } else {
                if (i >= len - clamp) return false;  // Normal: clamp is at end (3' end)
            }
            
            if (__builtin_expect((mismatches > max_mismatches), 0)) return false;
        }
    }
    
    return true;
}

// Comparison function for qsort
static int cmp_match_start(const void *a, const void *b) {
    const PrimerMatch *x = (const PrimerMatch*)a;
    const PrimerMatch *y = (const PrimerMatch*)b;
    return (x->start > y->start) - (x->start < y->start);
}

// Find all primer matches in a sequence (with SIMD optimization)
PrimerMatch* find_primer_matches(const char * restrict sequence, int seq_len,
                                 const uint8_t * restrict forward_mask, int f_len,
                                 const uint8_t * restrict reverse_mask, int r_len,
                                 const uint8_t * restrict forward_rc_mask,
                                 const uint8_t * restrict reverse_rc_mask,
                                 int max_mismatches, int clamp, int * restrict num_matches) {
    
    int capacity = 100;
    PrimerMatch *matches = (PrimerMatch*)malloc(capacity * sizeof(PrimerMatch));
    *num_matches = 0;
    
    // Build sequence mask once for entire sequence
    uint8_t *seq_mask_full = (uint8_t *)malloc(seq_len);
    build_seq_mask(sequence, seq_mask_full, seq_len);
    
    // Search for forward primer
    for (int pos = 0; pos <= seq_len - f_len; pos++) {
        if (matches_with_mismatches_simd(seq_mask_full + pos, forward_mask, f_len, 
                                         max_mismatches, clamp, false)) {
            if (__builtin_expect((*num_matches >= capacity), 0)) {
                capacity *= 2;
                matches = (PrimerMatch*)realloc(matches, capacity * sizeof(PrimerMatch));
            }
            
            matches[*num_matches].type = 0;
            matches[*num_matches].start = pos;
            matches[*num_matches].end = pos + f_len;
            strncpy(matches[*num_matches].primer_seq, sequence + pos, f_len);
            matches[*num_matches].primer_seq[f_len] = '\0';
            (*num_matches)++;
        }
    }
    
    // Search for forward primer RC
    for (int pos = 0; pos <= seq_len - f_len; pos++) {
        if (matches_with_mismatches_simd(seq_mask_full + pos, forward_rc_mask, f_len, 
                                         max_mismatches, clamp, true)) {
            if (__builtin_expect((*num_matches >= capacity), 0)) {
                capacity *= 2;
                matches = (PrimerMatch*)realloc(matches, capacity * sizeof(PrimerMatch));
            }
            
            matches[*num_matches].type = 1;
            matches[*num_matches].start = pos;
            matches[*num_matches].end = pos + f_len;
            strncpy(matches[*num_matches].primer_seq, sequence + pos, f_len);
            matches[*num_matches].primer_seq[f_len] = '\0';
            (*num_matches)++;
        }
    }
    
    // Search for reverse primer
    for (int pos = 0; pos <= seq_len - r_len; pos++) {
        if (matches_with_mismatches_simd(seq_mask_full + pos, reverse_mask, r_len, 
                                         max_mismatches, clamp, false)) {
            if (__builtin_expect((*num_matches >= capacity), 0)) {
                capacity *= 2;
                matches = (PrimerMatch*)realloc(matches, capacity * sizeof(PrimerMatch));
            }
            
            matches[*num_matches].type = 2;
            matches[*num_matches].start = pos;
            matches[*num_matches].end = pos + r_len;
            strncpy(matches[*num_matches].primer_seq, sequence + pos, r_len);
            matches[*num_matches].primer_seq[r_len] = '\0';
            (*num_matches)++;
        }
    }
    
    // Search for reverse primer RC
    for (int pos = 0; pos <= seq_len - r_len; pos++) {
        if (matches_with_mismatches_simd(seq_mask_full + pos, reverse_rc_mask, r_len, 
                                         max_mismatches, clamp, true)) {
            if (__builtin_expect((*num_matches >= capacity), 0)) {
                capacity *= 2;
                matches = (PrimerMatch*)realloc(matches, capacity * sizeof(PrimerMatch));
            }
            
            matches[*num_matches].type = 3;
            matches[*num_matches].start = pos;
            matches[*num_matches].end = pos + r_len;
            strncpy(matches[*num_matches].primer_seq, sequence + pos, r_len);
            matches[*num_matches].primer_seq[r_len] = '\0';
            (*num_matches)++;
        }
    }
    
    // Free the sequence mask
    free(seq_mask_full);
    
    // Sort matches by position using qsort
    if (*num_matches > 0) {
        qsort(matches, *num_matches, sizeof(PrimerMatch), cmp_match_start);
    }
    
    return matches;
}

// Free a single amplicon's memory
void free_amplicon(Amplicon *amp) {
    if (amp->sequence) free(amp->sequence);
    if (amp->source_file) free(amp->source_file);
    if (amp->seq_id) free(amp->seq_id);
    if (amp->forward_barcode) free(amp->forward_barcode);
    if (amp->reverse_barcode) free(amp->reverse_barcode);
}

// Write amplicons to temporary file and clear buffer
void write_amplicons_to_temp(AmpliconThreadData *data) {
    if (data->num_amplicons == 0) return;
    
    // Open temp file if not already open
    if (!data->temp_file.fp) {
        snprintf(data->temp_file.filepath, MAX_PATH_LENGTH, "%s/amplicons_t%d.tmp", 
                 data->temp_dir, data->thread_id);
        data->temp_file.fp = fopen(data->temp_file.filepath, "w");
        if (!data->temp_file.fp) {
            fprintf(stderr, "Error: Cannot create temporary file %s\n", data->temp_file.filepath);
            return;
        }
    }
    
    // Write amplicons to temp file (FASTA ONLY)
    for (int i = 0; i < data->num_amplicons; i++) {
        Amplicon *amp = &data->amplicons[i];
        
        fprintf(data->temp_file.fp, ">%s.source=%s.coordinates=%u-%u.orientation=%s",
               amp->seq_id, amp->source_file, amp->start, amp->end, amp->orientation);
        
        fprintf(data->temp_file.fp, ".fprimer=%s.rprimer=%s", 
                amp->forward_primer_seq, amp->reverse_primer_seq);
        
        if (amp->forward_barcode) {
            fprintf(data->temp_file.fp, ".fb=%s", amp->forward_barcode);
        }
        if (amp->reverse_barcode) {
            fprintf(data->temp_file.fp, ".rb=%s", amp->reverse_barcode);
        }
        fprintf(data->temp_file.fp, "\n");
        
        fprintf(data->temp_file.fp, "%s\n", amp->sequence);
        
        data->temp_file.count++;
    }
    
    // Flush to disk periodically
    if (data->temp_file.count % 1000 == 0) {
        fflush(data->temp_file.fp);
    }
    
    // Update total count
    data->total_amplicons_written += data->num_amplicons;
    
    // Free amplicon memory
    for (int i = 0; i < data->num_amplicons; i++) {
        free_amplicon(&data->amplicons[i]);
    }
    
    // Reset counter
    data->num_amplicons = 0;
}

// Process a sequence to find amplicons with aggressive memory management
void process_sequence_for_amplicons(const char * restrict seq_id, const char * restrict sequence, 
                                    const char * restrict source_file, 
                                    AmpliconThreadData * restrict data) {
    int seq_len = strlen(sequence);
    
    int num_matches;
    PrimerMatch *matches = find_primer_matches(sequence, seq_len,
                                               data->forward_mask, data->forward_len,
                                               data->reverse_mask, data->reverse_len,
                                               data->forward_rc_mask, data->reverse_rc_mask,
                                               data->mismatches, data->clamp, &num_matches);
    
    // Find valid amplicons
    for (int i = 0; i < num_matches - 1; i++) {
        if (matches[i].type != 0 && matches[i].type != 2) continue;
        
        for (int j = i + 1; j < num_matches; j++) {
            int length = matches[j].end - matches[i].start;
            
            if (length > data->max_length) break;
            
            // CRITICAL: Stop at next same-sense primer
            if (matches[j].type == 0 || matches[j].type == 2) break;
            
            if (length < data->min_length) continue;
            
            if (matches[j].type != 1 && matches[j].type != 3) continue;
            
            char orientation[3];
            if (matches[i].type == 0 && matches[j].type == 3) {
                strcpy(orientation, "FR");
            } else if (matches[i].type == 2 && matches[j].type == 1) {
                strcpy(orientation, "RF");
            } else if (matches[i].type == 0 && matches[j].type == 1) {
                strcpy(orientation, "FF");
            } else if (matches[i].type == 2 && matches[j].type == 3) {
                strcpy(orientation, "RR");
            } else {
                continue;
            }
            
            if (!data->include_offtarget && strcmp(orientation, "FR") != 0 && strcmp(orientation, "RF") != 0) {
                continue;
            }
            
            // Check if buffer is getting full and write immediately
            if (data->num_amplicons >= AMPLICON_BUFFER_SIZE) {
                write_amplicons_to_temp(data);
                // Force memory release after writing
#ifdef __linux__
                malloc_trim(0);
#endif
            }
            
            // Expand capacity if needed (shouldn't happen often with small buffer)
            if (__builtin_expect((data->num_amplicons >= data->capacity), 0)) {
                data->capacity += AMPLICON_BUFFER_SIZE;
                data->amplicons = (Amplicon*)realloc(data->amplicons, data->capacity * sizeof(Amplicon));
            }
            
            Amplicon *amp = &data->amplicons[data->num_amplicons];
            
            int amp_start = data->trim_primers ? matches[i].end : matches[i].start;
            int amp_end = data->trim_primers ? matches[j].start : matches[j].end;
            int amp_len = amp_end - amp_start;
            
            if (amp_len <= 0) {
                continue;
            }
            
            amp->sequence = (char*)malloc(amp_len + 1);
            strncpy(amp->sequence, sequence + amp_start, amp_len);
            amp->sequence[amp_len] = '\0';
            
            if (strcmp(orientation, "RF") == 0) {
                char *rc = reverse_complement_string(amp->sequence, amp_len);
                free(amp->sequence);
                amp->sequence = rc;
            }
            
            // Extract barcodes
            amp->forward_barcode = NULL;
            amp->reverse_barcode = NULL;
            
            if (strcmp(orientation, "FR") == 0) {
                if (data->fb_len > 0 && matches[i].start >= data->fb_len) {
                    amp->forward_barcode = (char*)malloc(data->fb_len + 1);
                    strncpy(amp->forward_barcode, sequence + matches[i].start - data->fb_len, data->fb_len);
                    amp->forward_barcode[data->fb_len] = '\0';
                }
                if (data->rb_len > 0 && matches[j].end + data->rb_len <= seq_len) {
                    amp->reverse_barcode = (char*)malloc(data->rb_len + 1);
                    strncpy(amp->reverse_barcode, sequence + matches[j].end, data->rb_len);
                    amp->reverse_barcode[data->rb_len] = '\0';
                }
            } else if (strcmp(orientation, "RF") == 0) {
                if (data->fb_len > 0 && matches[j].end + data->fb_len <= seq_len) {
                    char *temp = (char*)malloc(data->fb_len + 1);
                    strncpy(temp, sequence + matches[j].end, data->fb_len);
                    temp[data->fb_len] = '\0';
                    amp->forward_barcode = reverse_complement_string(temp, data->fb_len);
                    free(temp);
                }
                if (data->rb_len > 0 && matches[i].start >= data->rb_len) {
                    char *temp = (char*)malloc(data->rb_len + 1);
                    strncpy(temp, sequence + matches[i].start - data->rb_len, data->rb_len);
                    temp[data->rb_len] = '\0';
                    amp->reverse_barcode = reverse_complement_string(temp, data->rb_len);
                    free(temp);
                }
            }
            
            amp->source_file = strdup(source_file);
            amp->seq_id = strdup(seq_id);
            amp->start = matches[i].start;
            amp->end = matches[j].end;
            strcpy(amp->orientation, orientation);
            strcpy(amp->forward_primer_seq, matches[i].primer_seq);
            strcpy(amp->reverse_primer_seq, matches[j].primer_seq);
            
            data->num_amplicons++;
        }
    }
    
    free(matches);
}

// Worker thread for amplicon finding with TRUE STREAMING of 2bit files
void* amplicon_worker(void *arg) {
    AmpliconThreadData *data = (AmpliconThreadData*)arg;
    
    // Initialize temp file info
    data->temp_file.fp = NULL;
    data->temp_file.count = 0;
    data->total_amplicons_written = 0;
    
    // Process each file individually with streaming
    for (int idx = data->start_idx; idx < data->end_idx; idx++) {
        // Open file for streaming
        StreamingTwoBitFile seq_file = open_streaming_2bit(data->seq_files[idx]);
        if (!seq_file.fp) continue;
        
        char *base_name = basename(data->seq_files[idx]);
        
        // Process sequences one at a time
        TwoBitSequence seq;
        char *seq_buf = NULL;
        size_t seq_cap = 0;
        int seq_count = 0;
        
        while (read_next_2bit_sequence(&seq_file, &seq)) {
            // Allocate or resize buffer for this sequence
            if (seq.length + 1 > seq_cap) {
                seq_cap = seq.length + 1;
                seq_buf = (char*)realloc(seq_buf, seq_cap);
            }
            
            // Unpack and process
            unpack_sequence(seq.packed_data, seq.length, seq_buf);
            process_sequence_for_amplicons(seq.header, seq_buf, base_name, data);
            
            // Free the header for this sequence
            free_2bit_sequence(&seq);
            
            // Write amplicons frequently
            seq_count++;
            if (data->num_amplicons > 0 && (seq_count % 10 == 0)) {
                write_amplicons_to_temp(data);
                // Aggressive memory release every 10 sequences
#ifdef __linux__
                if (seq_count % 100 == 0) {
                    malloc_trim(0);
                }
#endif
            }
        }
        
        // Free sequence buffer
        if (seq_buf) {
            free(seq_buf);
            seq_buf = NULL;
        }
        
        // Write any remaining amplicons for this file
        if (data->num_amplicons > 0) {
            write_amplicons_to_temp(data);
        }
        
        // Close the streaming file
        close_streaming_2bit(&seq_file);
        
        // Force memory release after each file
#ifdef __linux__
        malloc_trim(0);
#endif
    }
    
    // Write any final remaining amplicons
    if (data->num_amplicons > 0) {
        write_amplicons_to_temp(data);
    }
    
    // Close temp file
    if (data->temp_file.fp) {
        fflush(data->temp_file.fp);
        fclose(data->temp_file.fp);
    }
    
    return NULL;
}

// Merge temporary files into final output
void merge_temp_files(AmpliconThreadData *thread_data, int num_threads, const char *output_file) {
    FILE *out = fopen(output_file, "w");
    if (!out) {
        fprintf(stderr, "Cannot create output file %s\n", output_file);
        return;
    }
    
    char buffer[MAX_LINE_LENGTH * 4];  // Buffer for efficient copying
    
    for (int t = 0; t < num_threads; t++) {
        if (thread_data[t].temp_file.count > 0) {
            FILE *temp = fopen(thread_data[t].temp_file.filepath, "r");
            if (temp) {
                while (fgets(buffer, sizeof(buffer), temp)) {
                    fputs(buffer, out);
                }
                fclose(temp);
                // Remove temporary file
                unlink(thread_data[t].temp_file.filepath);
            }
        }
    }
    
    fclose(out);
}

// Clean all temporary files created by this process
void cleanup_temp_files(const char *temp_dir) {
    DIR *dir = opendir(temp_dir);
    if (!dir) return;
    
    struct dirent *entry;
    char full_path[MAX_PATH_LENGTH];
    
    while ((entry = readdir(dir)) != NULL) {
        if (entry->d_name[0] == '.') continue;
        snprintf(full_path, sizeof(full_path), "%s/%s", temp_dir, entry->d_name);
        unlink(full_path);
    }
    
    closedir(dir);
    rmdir(temp_dir);
}

// Compress command with batching
int compress_command(const char *input_dir, const char *output_dir, int num_threads, size_t batch_size_mb) {
    // Create output directory
    mkdir(output_dir, 0755);
    
    // Find all sequence files
    int capacity = 1000;
    char **files = (char**)malloc(capacity * sizeof(char*));
    int num_files = 0;
    
    fprintf(stderr, "Scanning directory %s for FASTA files...\n", input_dir);
    find_sequence_files(input_dir, &files, &num_files, &capacity);
    fprintf(stderr, "Found %d FASTA files\n", num_files);
    
    if (num_files == 0) {
        fprintf(stderr, "No FASTA files found\n");
        free(files);
        return 1;
    }
    
    // Create temporary file for list
    char list_file[MAX_PATH_LENGTH];
    snprintf(list_file, sizeof(list_file), "%s/file_list.tmp", output_dir);
    FILE *list_fp = fopen(list_file, "w");
    if (!list_fp) {
        fprintf(stderr, "Cannot create file list\n");
        return 1;
    }
    
    pthread_mutex_t list_mutex = PTHREAD_MUTEX_INITIALIZER;
    
    // Create threads
    pthread_t *threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    CompressThreadData *thread_data = (CompressThreadData*)malloc(num_threads * sizeof(CompressThreadData));
    
    int files_per_thread = (num_files + num_threads - 1) / num_threads;
    size_t batch_size_bytes = batch_size_mb * 1024 * 1024;
    
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].thread_id = i;
        thread_data[i].files = files;
        thread_data[i].num_files = num_files;
        thread_data[i].start_idx = i * files_per_thread;
        thread_data[i].end_idx = ((i + 1) * files_per_thread < num_files) ? 
                                  (i + 1) * files_per_thread : num_files;
        thread_data[i].output_dir = (char*)output_dir;
        thread_data[i].list_fp = list_fp;
        thread_data[i].list_mutex = &list_mutex;
        thread_data[i].batch_size_bytes = batch_size_bytes;
        
        pthread_create(&threads[i], NULL, compress_worker, &thread_data[i]);
    }
    
    // Wait for threads
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }
    
    fclose(list_fp);
    
    // Output file list to stdout
    FILE *read_list = fopen(list_file, "r");
    if (read_list) {
        char line[MAX_PATH_LENGTH * 2];
        while (fgets(line, sizeof(line), read_list)) {
            printf("%s", line);
        }
        fclose(read_list);
        remove(list_file);
    }
    
    // Clean up
    for (int i = 0; i < num_files; i++) {
        free(files[i]);
    }
    free(files);
    free(threads);
    free(thread_data);
    pthread_mutex_destroy(&list_mutex);
    
    fprintf(stderr, "Compression complete. Created batched files (max %zu MB each)\n", batch_size_mb);
    fprintf(stderr, "Using IUPAC-aware encoding (degenerate bases mapped to lexicographically first valid base)\n");
    
    // Report memory usage
#ifdef __linux__
    FILE *status = fopen("/proc/self/status", "r");
    if (status) {
        char line[256];
        long vmpeak = 0;
        while (fgets(line, sizeof(line), status)) {
            if (sscanf(line, "VmPeak: %ld kB", &vmpeak) == 1) break;
        }
        fclose(status);
        if (vmpeak > 0) {
            fprintf(stderr, "Peak memory usage: %.2f MB\n", vmpeak / 1024.0);
        }
    }
#endif
    
    return 0;
}

// Run command with periodic amplicon writing
int run_command(const char *input_list, const char *primers_file, const char *output_file,
                int num_threads, int mismatches, int clamp, int min_length, int max_length,
                int fb_len, int rb_len, bool trim_primers, bool include_offtarget) {
    
    FILE *fp = fopen(primers_file, "r");
    if (!fp) {
        fprintf(stderr, "Cannot open primers file %s\n", primers_file);
        return 1;
    }
    
    char forward_primer[256], reverse_primer[256];
    if (fscanf(fp, "%s %s", forward_primer, reverse_primer) != 2) {
        fprintf(stderr, "Error reading primers from file\n");
        fclose(fp);
        return 1;
    }
    fclose(fp);
    
    fprintf(stderr, "Forward primer: %s\n", forward_primer);
    fprintf(stderr, "Reverse primer: %s\n", reverse_primer);
    
    int forward_len, reverse_len;
    uint8_t *forward_mask = primer_to_mask(forward_primer, &forward_len);
    uint8_t *reverse_mask = primer_to_mask(reverse_primer, &reverse_len);
    
    // Precompute RC masks once
    uint8_t *forward_rc_mask = reverse_complement_mask(forward_mask, forward_len);
    uint8_t *reverse_rc_mask = reverse_complement_mask(reverse_mask, reverse_len);
    
    int capacity = 1000;
    char **seq_files = (char**)malloc(capacity * sizeof(char*));
    int num_files = 0;
    
    fp = fopen(input_list, "r");
    if (!fp) {
        fprintf(stderr, "Cannot open input list %s\n", input_list);
        return 1;
    }
    
    char line[MAX_PATH_LENGTH * 2];
    while (fgets(line, sizeof(line), fp)) {
        if (num_files >= capacity) {
            capacity *= 2;
            seq_files = (char**)realloc(seq_files, capacity * sizeof(char*));
        }
        
        // Only read single column (no FASTQ support)
        line[strcspn(line, "\n\r\t")] = '\0';
        seq_files[num_files] = strdup(line);
        num_files++;
    }
    fclose(fp);
    
    fprintf(stderr, "Processing %d files with %d threads\n", num_files, num_threads);
    fprintf(stderr, "Using chunked reading mode (%d MB chunks, truly constant memory)\n", 
            MMAP_CHUNK_SIZE / (1024 * 1024));
    
    // Create temporary directory for intermediate files
    char temp_dir[MAX_PATH_LENGTH];
    snprintf(temp_dir, sizeof(temp_dir), "/tmp/amplicon_hunter_%d", getpid());
    mkdir(temp_dir, 0755);
    
    pthread_t *threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    AmpliconThreadData *thread_data = (AmpliconThreadData*)malloc(num_threads * sizeof(AmpliconThreadData));
    
    int files_per_thread = (num_files + num_threads - 1) / num_threads;
    
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].thread_id = i;
        thread_data[i].seq_files = seq_files;
        thread_data[i].num_files = num_files;
        thread_data[i].start_idx = i * files_per_thread;
        thread_data[i].end_idx = ((i + 1) * files_per_thread < num_files) ? 
                                  (i + 1) * files_per_thread : num_files;
        
        thread_data[i].forward_mask = forward_mask;
        thread_data[i].reverse_mask = reverse_mask;
        thread_data[i].forward_rc_mask = forward_rc_mask;
        thread_data[i].reverse_rc_mask = reverse_rc_mask;
        thread_data[i].forward_len = forward_len;
        thread_data[i].reverse_len = reverse_len;
        
        thread_data[i].mismatches = mismatches;
        thread_data[i].clamp = clamp;
        thread_data[i].min_length = min_length;
        thread_data[i].max_length = max_length;
        thread_data[i].fb_len = fb_len;
        thread_data[i].rb_len = rb_len;
        thread_data[i].trim_primers = trim_primers;
        thread_data[i].include_offtarget = include_offtarget;
        
        // Initialize amplicon buffer with SMALL initial capacity for streaming
        thread_data[i].amplicons = (Amplicon*)malloc(AMPLICON_BUFFER_SIZE * sizeof(Amplicon));
        thread_data[i].num_amplicons = 0;
        thread_data[i].capacity = AMPLICON_BUFFER_SIZE;
        thread_data[i].temp_dir = temp_dir;
        
        pthread_create(&threads[i], NULL, amplicon_worker, &thread_data[i]);
    }
    
    // Wait for threads
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }
    
    // Count total amplicons
    int total_amplicons = 0;
    for (int i = 0; i < num_threads; i++) {
        total_amplicons += thread_data[i].total_amplicons_written;
    }
    
    fprintf(stderr, "Found %d amplicons\n", total_amplicons);
    
    // Get memory usage stats if available
#ifdef __linux__
    FILE *status = fopen("/proc/self/status", "r");
    if (status) {
        char line[256];
        long vmrss = 0, vmpeak = 0;
        while (fgets(line, sizeof(line), status)) {
            if (sscanf(line, "VmRSS: %ld kB", &vmrss) == 1) continue;
            if (sscanf(line, "VmPeak: %ld kB", &vmpeak) == 1) continue;
        }
        fclose(status);
        if (vmpeak > 0) {
            fprintf(stderr, "Peak memory usage: %.2f MB\n", vmpeak / 1024.0);
        }
    }
#endif
    
    // Merge temporary files into final output
    merge_temp_files(thread_data, num_threads, output_file);
    
    // Clean up temporary directory
    cleanup_temp_files(temp_dir);
    
    // Clean up
    for (int i = 0; i < num_threads; i++) {
        free(thread_data[i].amplicons);
    }
    
    free(forward_mask);
    free(reverse_mask);
    free(forward_rc_mask);
    free(reverse_rc_mask);
    
    for (int i = 0; i < num_files; i++) {
        free(seq_files[i]);
    }
    free(seq_files);
    free(threads);
    free(thread_data);
    
    return 0;
}

void print_usage() {
    printf("AmpliconHunter v%s - IUPAC-aware, chunked-reading amplicon finder\n\n", VERSION);
    printf("Uses %d MB chunks for truly constant memory usage\n\n", MMAP_CHUNK_SIZE / (1024 * 1024));
    printf("Usage:\n");
    printf("  amplicon_hunter compress --input-dir <dir> --output <dir> [--threads <n>] [--batch-size <MB>]\n");
    printf("  amplicon_hunter run --input <file_list> --primers <file> --output <file> [options]\n\n");
    
    printf("Compress options:\n");
    printf("  --input-dir <dir>    Directory containing FASTA files\n");
    printf("  --output <dir>       Output directory for compressed files\n");
    printf("  --threads <n>        Number of threads (default: 8)\n");
    printf("  --batch-size <MB>    Max size per batch file in MB (default: 500)\n\n");
    
    printf("Run options:\n");
    printf("  --input <file>       File list from compress (single column)\n");
    printf("  --primers <file>     Primers file (forward reverse)\n");
    printf("  --output <file>      Output FASTA file\n");
    printf("  --threads <n>        Number of threads (default: 8)\n");
    printf("  --mismatches <n>     Maximum mismatches (default: 0)\n");
    printf("  --clamp <n>          3' clamp size (default: 5)\n");
    printf("  --min-length <n>     Minimum amplicon length (default: 50)\n");
    printf("  --max-length <n>     Maximum amplicon length (default: 5000)\n");
    printf("  --fb-len <n>         Forward barcode length (default: 0)\n");
    printf("  --rb-len <n>         Reverse barcode length (default: 0)\n");
    printf("  --trim-primers       Trim primer sequences\n");
    printf("  --include-offtarget  Include FF and RR orientations\n");
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        print_usage();
        return 1;
    }
    
    init_lookup_tables();
    
    if (strcmp(argv[1], "compress") == 0) {
        char *input_dir = NULL;
        char *output_dir = NULL;
        int num_threads = DEFAULT_THREADS;
        size_t batch_size_mb = DEFAULT_BATCH_SIZE_MB;
        
        for (int i = 2; i < argc; i++) {
            if (strcmp(argv[i], "--input-dir") == 0 && i + 1 < argc) {
                input_dir = argv[++i];
            } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
                output_dir = argv[++i];
            } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
                num_threads = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--batch-size") == 0 && i + 1 < argc) {
                batch_size_mb = atoi(argv[++i]);
            }
        }
        
        if (!input_dir || !output_dir) {
            fprintf(stderr, "Error: Missing required arguments\n");
            print_usage();
            return 1;
        }
        
        return compress_command(input_dir, output_dir, num_threads, batch_size_mb);
        
    } else if (strcmp(argv[1], "run") == 0) {
        char *input_list = NULL;
        char *primers_file = NULL;
        char *output_file = NULL;
        int num_threads = DEFAULT_THREADS;
        int mismatches = 0;
        int clamp = 5;
        int min_length = 50;
        int max_length = 5000;
        int fb_len = 0;
        int rb_len = 0;
        bool trim_primers = false;
        bool include_offtarget = false;
        
        for (int i = 2; i < argc; i++) {
            if (strcmp(argv[i], "--input") == 0 && i + 1 < argc) {
                input_list = argv[++i];
            } else if (strcmp(argv[i], "--primers") == 0 && i + 1 < argc) {
                primers_file = argv[++i];
            } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
                output_file = argv[++i];
            } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
                num_threads = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--mismatches") == 0 && i + 1 < argc) {
                mismatches = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--clamp") == 0 && i + 1 < argc) {
                clamp = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--min-length") == 0 && i + 1 < argc) {
                min_length = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--max-length") == 0 && i + 1 < argc) {
                max_length = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--fb-len") == 0 && i + 1 < argc) {
                fb_len = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--rb-len") == 0 && i + 1 < argc) {
                rb_len = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--trim-primers") == 0) {
                trim_primers = true;
            } else if (strcmp(argv[i], "--include-offtarget") == 0) {
                include_offtarget = true;
            }
        }
        
        if (!input_list || !primers_file || !output_file) {
            fprintf(stderr, "Error: Missing required arguments\n");
            print_usage();
            return 1;
        }
        
        return run_command(input_list, primers_file, output_file, num_threads,
                          mismatches, clamp, min_length, max_length,
                          fb_len, rb_len, trim_primers, include_offtarget);
        
    } else {
        fprintf(stderr, "Error: Unknown command '%s'\n", argv[1]);
        print_usage();
        return 1;
    }
    
    return 0;
}