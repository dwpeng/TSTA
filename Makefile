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
SRCS := $(wildcard src/tsta_*.c)
OBJS := $(SRCS:.c=.o)
SHARED_OBJS := $(SRCS:.c=.shared.o)

.PHONY: all clean shared static install uninstall check test

all: static shared

static: $(STATIC_LIB)

shared: $(SHARED_LIB_FULL)

$(STATIC_LIB): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

$(SHARED_LIB_FULL): $(SHARED_OBJS)
	$(CC) -shared -Wl,-soname,$(SHARED_LIB_MAJOR) -o $@ $^ -lpthread
	ln -sf $(SHARED_LIB_FULL) $(SHARED_LIB_MAJOR)
	ln -sf $(SHARED_LIB_MAJOR) $(SHARED_LIB)

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

src/%.shared.o: src/%.c
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

clean:
	rm -f src/tsta_*.o src/tsta_*.shared.o $(STATIC_LIB) $(SHARED_LIB) $(SHARED_LIB_MAJOR) $(SHARED_LIB_FULL) test/test_main

test: $(STATIC_LIB)
	$(CC) $(CFLAGS) -o test/test_main test/main.c $(STATIC_LIB)

check: test
	./test/test_main

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
