CC = gcc
CFLAGS = -Wall -O3 -march=native -mavx2 -mpopcnt -pthread -ffast-math -flto -fno-math-errno -fno-trapping-math -DNDEBUG
LDFLAGS = -pthread -lm -flto
TARGET = amplicon_hunter

# Source files
SRCS = amplicon_hunter.c
OBJS = $(SRCS:.c=.o)

# Default target - optimized for streaming memory
all: $(TARGET)

# Build the executable
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Debug build with address sanitizer for memory leak detection
debug: CFLAGS = -Wall -g -O0 -pthread -DDEBUG -fsanitize=address -fno-omit-frame-pointer
debug: LDFLAGS = -pthread -lm -fsanitize=address
debug: clean $(TARGET)

# Memory-optimized build with aggressive optimization
memory-opt: CFLAGS = -Wall -O3 -march=native -mavx2 -mpopcnt -pthread -ffast-math -flto -fno-math-errno -fno-trapping-math -DNDEBUG -fdata-sections -ffunction-sections
memory-opt: LDFLAGS = -pthread -lm -flto -Wl,--gc-sections
memory-opt: clean $(TARGET)
	strip $(TARGET)

# Debug build with thread sanitizer
debug-thread: CFLAGS = -Wall -g -O0 -pthread -DDEBUG -fsanitize=thread
debug-thread: LDFLAGS = -pthread -lm -fsanitize=thread
debug-thread: clean $(TARGET)

# Profile-guided optimization
pgo: clean
	$(CC) $(CFLAGS) -fprofile-generate -o $(TARGET) $(SRCS) $(LDFLAGS)
	@echo "Run typical workloads, then 'make pgo-use'"

pgo-use:
	$(CC) $(CFLAGS) -fprofile-use -o $(TARGET) $(SRCS) $(LDFLAGS)

# Ultra-optimized build for maximum performance
ultra: CFLAGS = -Wall -Ofast -march=native -mavx2 -mpopcnt -pthread -ffast-math -funroll-loops -flto -fno-math-errno -fno-trapping-math -DNDEBUG
ultra: LDFLAGS = -pthread -lm -flto
ultra: clean $(TARGET)

# Clean all temp files and build artifacts
clean:
	rm -f $(TARGET) $(OBJS) *.gcda *.gcno
	rm -rf /tmp/amplicon_hunter_*
	@echo "Cleaned all temporary files and build artifacts"

# Deep clean - also removes any leftover batch files
deep-clean: clean
	find . -name "batch_t*.2bit" -delete
	find . -name "batch_t*.qual" -delete
	find . -name "*.tmp" -delete
	@echo "Deep clean completed"

# Install
install: $(TARGET)
	install -m 755 $(TARGET) /usr/local/bin/

# Uninstall
uninstall:
	rm -f /usr/local/bin/$(TARGET)

# Test compress (FASTA only)
test-compress:
	./$(TARGET) compress --input-dir test/fasta --output test/compressed --threads 8 --batch-size 500

# Test run
test-run:
	./$(TARGET) run --input test/file_list.txt --primers test/primers.txt \
		--output test_output.fa --threads 8 --mismatches 2 --clamp 3 \
		--min-length 1000 --max-length 2000

# Test with large dataset (stress test for memory)
test-large:
	./$(TARGET) run --input test/large_file_list.txt --primers test/primers.txt \
		--output test_large_output.fa --threads 16 --mismatches 3 --clamp 3 \
		--min-length 50 --max-length 10000 --include-offtarget

# Memory test with valgrind
memtest: debug
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes \
		./$(TARGET) run --input test/file_list.txt --primers test/primers.txt \
		--output test_memory.fa --threads 4 --mismatches 2 --clamp 3

# Memory profiling with massif
massif: memory-opt
	valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out \
		./$(TARGET) run --input test/file_list.txt --primers test/primers.txt \
		--output test_massif.fa --threads 4 --mismatches 2 --clamp 3
	ms_print massif.out > massif_report.txt
	@echo "Memory profile saved to massif_report.txt"

# Benchmark with different thread counts
benchmark: ultra
	@echo "Running benchmark with different thread counts..."
	@for t in 1 2 4 8 16; do \
		echo "Threads: $$t"; \
		rm -rf /tmp/amplicon_hunter_*; \
		time ./$(TARGET) run --input test/file_list.txt --primers test/primers.txt \
			--output /dev/null --threads $$t --mismatches 2 --clamp 3 \
			--min-length 1000 --max-length 2000 2>&1 | grep real; \
	done

# Performance profiling with perf
profile: ultra
	perf record -g ./$(TARGET) run --input test/file_list.txt --primers test/primers.txt \
		--output test_profile.fa --threads 8 --mismatches 2 --clamp 3
	perf report

# Check for memory leaks in compress
test-compress-memory: debug
	valgrind --leak-check=full ./$(TARGET) compress --input-dir test/fasta \
		--output test/compressed_memtest --threads 4 --batch-size 100

# Help
help:
	@echo "Available targets:"
	@echo "  all              - Build the standard optimized version"
	@echo "  memory-opt       - Build with aggressive memory optimization"
	@echo "  debug            - Build with debug symbols and AddressSanitizer"
	@echo "  debug-thread     - Build with debug symbols and ThreadSanitizer"
	@echo "  ultra            - Build with maximum performance optimizations"
	@echo "  pgo              - Build for profile-guided optimization"
	@echo "  clean            - Remove build artifacts and temp files"
	@echo "  deep-clean       - Remove all artifacts including batch files"
	@echo "  install          - Install to /usr/local/bin"
	@echo "  test-compress    - Test compression functionality"
	@echo "  test-run         - Test amplicon finding"
	@echo "  memtest          - Run with valgrind to check for memory leaks"
	@echo "  massif           - Profile memory usage with valgrind massif"
	@echo "  benchmark        - Run performance benchmarks"
	@echo "  profile          - Run with perf profiling"

.PHONY: all clean deep-clean debug debug-thread memory-opt ultra pgo pgo-use \
        install uninstall test-compress test-run test-large memtest massif \
        benchmark profile help test-compress-memory