CC ?= gcc
AR ?= ar
ARFLAGS ?= rcs
CFLAGS ?= -O3 -g -msse4.2 -I.
ZLIB_LIBS ?=

# Library version
VERSION_MAJOR := 1
VERSION_MINOR := 0
VERSION_PATCH := 0
VERSION := $(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)

ZLIB_CHECK := $(shell printf '#include <zlib.h>\nint main(void){return 0;}\n' | $(CC) -x c - -lz -o /tmp/tsta_zlib_check >/dev/null 2>&1 && echo yes)
ifeq ($(ZLIB_CHECK),yes)
	CFLAGS += -DTSTA_HAVE_ZLIB
	ZLIB_LIBS += -lz
endif

# Library names
STATIC_LIB := libtsta.a
SHARED_LIB := libtsta.so
SHARED_LIB_MAJOR := $(SHARED_LIB).$(VERSION_MAJOR)
SHARED_LIB_FULL := $(SHARED_LIB_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)

INCLUDE_DIRS := -I./include
CFLAGS += $(INCLUDE_DIRS)
SRCS := $(wildcard src/tsta_*.c)
OBJS := $(SRCS:.c=.o)
SHARED_OBJS := $(SRCS:.c=.shared.o)

.PHONY: all clean shared static install uninstall

all: static shared

static: $(STATIC_LIB)

shared: $(SHARED_LIB_FULL)

$(STATIC_LIB): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

$(SHARED_LIB_FULL): $(SHARED_OBJS)
	$(CC) -shared -Wl,-soname,$(SHARED_LIB_MAJOR) -o $@ $^ $(ZLIB_LIBS) -lpthread
	ln -sf $(SHARED_LIB_FULL) $(SHARED_LIB_MAJOR)
	ln -sf $(SHARED_LIB_MAJOR) $(SHARED_LIB)

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

src/%.shared.o: src/%.c
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

clean:
	rm -f src/tsta_*.o src/tsta_*.shared.o $(STATIC_LIB) $(SHARED_LIB) $(SHARED_LIB_MAJOR) $(SHARED_LIB_FULL) test/test_main

test: $(STATIC_LIB)
	$(CC) $(CFLAGS) -o test/test_main test/main.c $(STATIC_LIB) $(ZLIB_LIBS)

# Installation
PREFIX ?= /usr/local
INCLUDE_DIR ?= $(PREFIX)/include
LIB_DIR ?= $(PREFIX)/lib

install: all
	install -d $(DESTDIR)$(INCLUDE_DIR)/tsta
	install -m 644 include/tsta.h $(DESTDIR)$(INCLUDE_DIR)/
	install -m 644 include/tsta_api.h $(DESTDIR)$(INCLUDE_DIR)/tsta/
	install -m 644 $(STATIC_LIB) $(DESTDIR)$(LIB_DIR)/
	install -m 755 $(SHARED_LIB_FULL) $(DESTDIR)$(LIB_DIR)/
	cd $(DESTDIR)$(LIB_DIR) && ln -sf $(SHARED_LIB_FULL) $(SHARED_LIB_MAJOR)
	cd $(DESTDIR)$(LIB_DIR) && ln -sf $(SHARED_LIB_MAJOR) $(SHARED_LIB)
	ldconfig || true

install-static: static
	install -d $(DESTDIR)$(INCLUDE_DIR)/tsta
	install -m 644 include/tsta.h $(DESTDIR)$(INCLUDE_DIR)/
	install -m 644 include/tsta_api.h $(DESTDIR)$(INCLUDE_DIR)/tsta/
	install -m 644 $(STATIC_LIB) $(DESTDIR)$(LIB_DIR)/

install-shared: shared
	install -d $(DESTDIR)$(INCLUDE_DIR)/tsta
	install -m 644 include/tsta.h $(DESTDIR)$(INCLUDE_DIR)/
	install -m 644 include/tsta_api.h $(DESTDIR)$(INCLUDE_DIR)/tsta/
	install -m 755 $(SHARED_LIB_FULL) $(DESTDIR)$(LIB_DIR)/
	cd $(DESTDIR)$(LIB_DIR) && ln -sf $(SHARED_LIB_FULL) $(SHARED_LIB_MAJOR)
	cd $(DESTDIR)$(LIB_DIR) && ln -sf $(SHARED_LIB_MAJOR) $(SHARED_LIB)
	ldconfig || true

uninstall:
	rm -f $(DESTDIR)$(INCLUDE_DIR)/tsta.h
	rm -rf $(DESTDIR)$(INCLUDE_DIR)/tsta
	rm -f $(DESTDIR)$(LIB_DIR)/$(STATIC_LIB)
	rm -f $(DESTDIR)$(LIB_DIR)/$(SHARED_LIB)
	rm -f $(DESTDIR)$(LIB_DIR)/$(SHARED_LIB_MAJOR)
	rm -f $(DESTDIR)$(LIB_DIR)/$(SHARED_LIB_FULL)
	ldconfig || true