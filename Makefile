CC ?= gcc
AR ?= ar
ARFLAGS ?= rcs
CFLAGS ?= -O3 -g -march=native -fvisibility=hidden

# Library version
VERSION_MAJOR := 1
VERSION_MINOR := 0
VERSION_PATCH := 0
VERSION := $(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)

# Library names
STATIC_LIB := libtsta.a
SHARED_LIB := libtsta.so
SHARED_LIB_MAJOR := $(SHARED_LIB).$(VERSION_MAJOR)
SHARED_LIB_FULL := $(SHARED_LIB_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)

INCLUDE_DIRS := -I./src -I./include
CFLAGS += $(INCLUDE_DIRS)
# Keep the standalone SIMD-only module out of the threaded main library.
SRCS := $(filter-out src/tsta_psa_simd.c,$(wildcard src/tsta_*.c))
OBJS := $(SRCS:.c=.o)
SHARED_OBJS := $(SRCS:.c=.shared.o)

# ── Sanitizer support ──────────────────────────────────────────────────
# Usage:
#   make clean all SANITIZE=asan       # address + leak
#   make clean all SANITIZE=ubsan      # undefined behavior
#   make clean all SANITIZE=all        # both
#   make check_sanitize                # rebuild + test with asan+ubsan
#
# Also respected: ASAN_OPTIONS, UBSAN_OPTIONS (passed to test runs).

SANITIZE ?=

ifeq ($(SANITIZE),1)
  SANITIZE_FLAGS := -fsanitize=address,leak,undefined
  SANITIZE_ASAN := 1
  SANITIZE_UBSAN := 1
else ifeq ($(SANITIZE),asan)
  SANITIZE_FLAGS := -fsanitize=address,leak
  SANITIZE_ASAN := 1
else ifeq ($(SANITIZE),ubsan)
  SANITIZE_FLAGS := -fsanitize=undefined
  SANITIZE_UBSAN := 1
else ifeq ($(SANITIZE),all)
  SANITIZE_FLAGS := -fsanitize=address,leak,undefined
  SANITIZE_ASAN := 1
  SANITIZE_UBSAN := 1
endif

ifneq ($(SANITIZE),)
  SANITIZE_FLAGS += -fno-omit-frame-pointer -fno-sanitize-recover=all
  override CFLAGS := -O0 -g $(SANITIZE_FLAGS) $(filter-out -O3 -O2 -O1 -g,$(CFLAGS))
  override CXXFLAGS := -O0 -g $(SANITIZE_FLAGS) $(filter-out -O3 -O2 -O1 -g,$(CXXFLAGS))

  # Test runner for sanitized builds
  define run_sanitized
    $(if $(SANITIZE_ASAN),ASAN_OPTIONS=detect_leaks=1:abort_on_error=1)
    $(if $(SANITIZE_UBSAN),UBSAN_OPTIONS=halt_on_error=1)
    $(1)
  endef
else
  TEST_RUNNER :=
endif

.PHONY: all clean shared static install uninstall check test check_cpp examples \
        check_examples check_sanitize psa_simd test_psa_simd check_psa_simd \
        verify_psa_simd_nopthread format check-format install-hooks

all: static shared

static: $(STATIC_LIB)

shared: $(SHARED_LIB_FULL)

$(STATIC_LIB): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

$(SHARED_LIB_FULL): $(SHARED_OBJS)
	$(CC) -shared -Wl,-soname,$(SHARED_LIB_MAJOR) -o $@ $^ -lpthread
	ln -sf $(SHARED_LIB_FULL) $(SHARED_LIB_MAJOR)
	ln -sf $(SHARED_LIB_MAJOR) $(SHARED_LIB)

# ── Standalone SIMD-only PSA module (single-threaded, no pthread) ───────
# Builds and links with no -lpthread. Output is bit-identical to the
# threaded PSA for the same configuration and inputs.

PSA_SIMD_STATIC_LIB   := libtsta_psa_simd.a
PSA_SIMD_SHARED_MAJOR := libtsta_psa_simd.so.1
PSA_SIMD_SHARED_FULL  := $(PSA_SIMD_SHARED_MAJOR).0.0

# The module reuses the common helpers (config/state/result/aligned alloc),
# so the standalone archive bundles tsta_common.o alongside it.
$(PSA_SIMD_STATIC_LIB): src/tsta_common.o src/tsta_psa_simd.o
	$(AR) $(ARFLAGS) $@ $^

$(PSA_SIMD_SHARED_FULL): src/tsta_common.shared.o src/tsta_psa_simd.shared.o
	$(CC) -shared -Wl,-soname,$(PSA_SIMD_SHARED_MAJOR) -o $@ $^
	ln -sf $(PSA_SIMD_SHARED_FULL) $(PSA_SIMD_SHARED_MAJOR)
	ln -sf $(PSA_SIMD_SHARED_MAJOR) libtsta_psa_simd.so

psa_simd: $(PSA_SIMD_STATIC_LIB) $(PSA_SIMD_SHARED_FULL)

test_psa_simd: $(STATIC_LIB) $(PSA_SIMD_STATIC_LIB)
	$(CC) $(CFLAGS) -o test/test_psa_simd test/test_psa_simd.c $(STATIC_LIB) $(PSA_SIMD_STATIC_LIB)

check_psa_simd: test_psa_simd
	$(if $(SANITIZE),$(call run_sanitized,./test/test_psa_simd),./test/test_psa_simd)

verify_psa_simd_nopthread: $(PSA_SIMD_STATIC_LIB) $(PSA_SIMD_SHARED_FULL)
	@if nm $(PSA_SIMD_STATIC_LIB) | grep -Ei 'pthread|threadpool'; then \
	  echo "FAIL: pthread/threadpool symbols in $(PSA_SIMD_STATIC_LIB)"; exit 1; \
	else \
	  echo "OK: no pthread/threadpool symbols in $(PSA_SIMD_STATIC_LIB)"; fi
	@if nm -D $(PSA_SIMD_SHARED_FULL) | grep -Ei 'pthread|threadpool'; then \
	  echo "FAIL: pthread/threadpool symbols in $(PSA_SIMD_SHARED_FULL)"; exit 1; \
	else \
	  echo "OK: no pthread/threadpool symbols in $(PSA_SIMD_SHARED_FULL)"; fi

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

src/%.shared.o: src/%.c
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

CXX ?= g++
CXXFLAGS ?= -O3 -g -march=native -fvisibility=hidden

test_cpp: $(STATIC_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -std=c++17 -o test/test_cpp test/main.cpp $(STATIC_LIB)

check_cpp: test_cpp
	$(if $(SANITIZE),$(call run_sanitized,./test/test_cpp),./test/test_cpp)

# ── Formatting (clang-format) ──────────────────────────────────────────

FORMAT_FILES := $(wildcard src/*.c src/*.h include/*.h test/*.c test/*.cpp)

format:
	clang-format -i $(FORMAT_FILES)

check-format:
	clang-format --dry-run --Werror $(FORMAT_FILES)

# Use the project-local hooks (see .githooks/pre-commit).
install-hooks:
	git config core.hooksPath .githooks

clean:
	rm -f src/tsta_*.o src/tsta_*.shared.o $(STATIC_LIB) $(SHARED_LIB) $(SHARED_LIB_MAJOR) $(SHARED_LIB_FULL) test/test_main test/test_cpp examples/example_c examples/example_cpp $(PSA_SIMD_STATIC_LIB) libtsta_psa_simd.so $(PSA_SIMD_SHARED_MAJOR) $(PSA_SIMD_SHARED_FULL) test/test_psa_simd

test: $(STATIC_LIB)
	$(CC) $(CFLAGS) -o test/test_main test/main.c $(STATIC_LIB)

check: test
	$(if $(SANITIZE),$(call run_sanitized,./test/test_main),./test/test_main)

check_sanitize:
	@$(MAKE) clean
	@$(MAKE) check SANITIZE=all
	@$(MAKE) check_cpp SANITIZE=all
	@$(MAKE) check_examples SANITIZE=all
	@echo "=== All sanitizer checks passed ==="

# ── Examples ─────────────────────────────────────────────────────────

examples: $(STATIC_LIB)
	$(CC) $(CFLAGS) -o examples/example_c examples/main.c $(STATIC_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -std=c++17 -o examples/example_cpp examples/main.cc $(STATIC_LIB)

check_examples: examples
	$(if $(SANITIZE),$(call run_sanitized,./examples/example_c),./examples/example_c)
	$(if $(SANITIZE),$(call run_sanitized,./examples/example_cpp),./examples/example_cpp)

# Installation
PREFIX ?= /usr/local
INCLUDE_DIR ?= $(PREFIX)/include
LIB_DIR ?= $(PREFIX)/lib

install: all
	install -d $(DESTDIR)$(INCLUDE_DIR)
	install -m 644 include/tsta.h $(DESTDIR)$(INCLUDE_DIR)/
	install -m 644 $(STATIC_LIB) $(DESTDIR)$(LIB_DIR)/
	install -m 755 $(SHARED_LIB_FULL) $(DESTDIR)$(LIB_DIR)/
	cd $(DESTDIR)$(LIB_DIR) && ln -sf $(SHARED_LIB_FULL) $(SHARED_LIB_MAJOR)
	cd $(DESTDIR)$(LIB_DIR) && ln -sf $(SHARED_LIB_MAJOR) $(SHARED_LIB)
	ldconfig || true

install-static: static
	install -d $(DESTDIR)$(INCLUDE_DIR)
	install -m 644 include/tsta.h $(DESTDIR)$(INCLUDE_DIR)/
	install -m 644 $(STATIC_LIB) $(DESTDIR)$(LIB_DIR)/

install-shared: shared
	install -d $(DESTDIR)$(INCLUDE_DIR)
	install -m 644 include/tsta.h $(DESTDIR)$(INCLUDE_DIR)/
	install -m 755 $(SHARED_LIB_FULL) $(DESTDIR)$(LIB_DIR)/
	cd $(DESTDIR)$(LIB_DIR) && ln -sf $(SHARED_LIB_FULL) $(SHARED_LIB_MAJOR)
	cd $(DESTDIR)$(LIB_DIR) && ln -sf $(SHARED_LIB_MAJOR) $(SHARED_LIB)
	ldconfig || true

uninstall:
	rm -f $(DESTDIR)$(INCLUDE_DIR)/tsta.h
	rm -f $(DESTDIR)$(LIB_DIR)/$(STATIC_LIB)
	rm -f $(DESTDIR)$(LIB_DIR)/$(SHARED_LIB)
	rm -f $(DESTDIR)$(LIB_DIR)/$(SHARED_LIB_MAJOR)
	rm -f $(DESTDIR)$(LIB_DIR)/$(SHARED_LIB_FULL)
	ldconfig || true
