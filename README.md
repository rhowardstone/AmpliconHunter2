# AmpliconHunter2

SIMD-accelerated, in-silico PCR for very large FASTA sets. AVX2 bit-mask IUPAC matching, constant-memory streaming with 2-bit batches, multi-threaded amplicon calling with optional primer trimming and barcode extraction.

See our [webserver](https://ah2.engr.uconn.edu/) for an easy-to-use visual interface running AmpliconHunter2 on five genomes DBs of various sizes (~5.7k genomes to ~2.4M).

## Features

* AVX2 one-byte IUPAC masks with exact 3′ clamp and up to `--mismatches` substitutions.
* Streaming 2-bit batches and mmap-aware reads for low RSS on huge inputs.
* FR and RF by default. `--include-offtarget` adds FF and RR. RF is reverse-complemented on output.
* Melting temperature annotation and optional filtering (`--Tm`).
* Optional fixed-length flank barcodes (`--fb-len`, `--rb-len`) and `--trim-primers`.
* FASTA-only; see [AmpliconHunter Version 1.1](https://github.com/rhowardstone/AmpliconHunter) if you require more advanced features.

## Build

```bash
git clone https://github.com/rhowardstone/AmpliconHunter2.git
cd AmpliconHunter2; make;
sudo cp amplicon_hunter /usr/local/bin/  #for local installation
```

CPU must support AVX2. GCC 9+ recommended.

## Quick start

1. Batch-compress FASTA to 2-bit and capture the file list for the test data:

```bash
amplicon_hunter compress --input-dir test/fasta --output test/compressed
ls -d "$PWD"/test/compressed/*.2bit > test/file_list.txt
```


2. Prepare `primers.txt` on two lines: `<FORWARD>\n<REVERSE>` (5' -> 3' direction), test/V1V9_primers.txt contains the following for V1V9 (full-length 16S):

```
AGRGTTYGATYMTGGCTCAG
RGYTACCTTGTTACGACTT
```

3. Run the finder:

```bash
amplicon_hunter run \
  --input test/file_list.txt \
  --primers test/V1V9_primers.txt \
  --output test/amplicons.fa

diff test/amplicons{,-check}.fa | wc -l  #should show zero
```

## Output

FASTA with annotated headers:

```
>seqid.fileID=original_file.fa.source=batch_t0_b1.2bit.coordinates=12345-13567.Tm=60.54.orientation=FR.fprimer=... .rprimer=... .fb=ACGT.rb=TGCA
ACGT...
```

Header fields: source batch file, genomic coordinates, melting temperature, orientation, matched primer snippets, and optional barcodes.


## Full syntax:
```
Usage:
  amplicon_hunter compress --input-dir <dir> --output <dir> [--threads <n>] [--batch-size <MB>]
  amplicon_hunter run --input <file_list> --primers <file> --output <file> [options]

Compress options:
  --input-dir <dir>    Directory containing FASTA files
  --output <dir>       Output directory for compressed files
  --threads <n>        Number of threads (default: 8)
  --batch-size <MB>    Max size per batch file in MB (default: 500)

Run options:
  --input <file>       File list from compress (single column)
  --primers <file>     Primers file (forward reverse)
  --output <file>      Output FASTA file
  --threads <n>        Number of threads (default: 8)
  --mismatches <n>     Maximum mismatches (default: 0)
  --clamp <n>          3' clamp size (default: 5)
  --min-length <n>     Minimum amplicon length (default: 50)
  --max-length <n>     Maximum amplicon length (default: 5000)
  --fb-len <n>         Forward barcode length (default: 0)
  --rb-len <n>         Reverse barcode length (default: 0)
  --trim-primers       Trim primer sequences
  --include-offtarget  Include FF and RR orientations

Tm calculation options:
  --Tm <n>             Minimum melting temperature threshold (default: 0, no filtering)
  --dnac1 <n>          Primer concentration in nM (default: 1000)
  --dnac2 <n>          Template concentration in nM (default: 25)
  --Na <n>             Sodium concentration in mM (default: 50)
  --K <n>              Potassium concentration in mM (default: 0)
  --Tris <n>           Tris buffer concentration in mM (default: 0)
  --Mg <n>             Magnesium concentration in mM (default: 4)
  --dNTPs <n>          dNTP concentration in mM (default: 1.6)
  --saltcorr <n>       Salt correction method 0-7 (default: 5)
```


## Notes and limits

* Mismatches are substitutions only; No indels.
* Non-ACGT input symbols are mapped to the alphabetically first possible base during 2-bit packing. Uppercase and lowercase treated identically.
* Temp files are written under `/tmp/amplicon_hunter_<pid>` and merged automatically.
* If your input is FASTQ format, please see the Python version of our tool: [AmpliconHunter Version 1.1](https://github.com/rhowardstone/AmpliconHunter)

## Benchmarking repo

For detailed usage examples, and for comparing our intermediate versions, please see our [Benchmarking repo for AmpliconHunter2](https://github.com/rhowardstone/AmpliconHunter2_benchmark).

## Cite

If this tool supports your work, cite the AmpliconHunter2 application note (will update with DOI when available).

## License
AmpliconHunter2 is licensed under the MIT License. See LICENSE for details.
