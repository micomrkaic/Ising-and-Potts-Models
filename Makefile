# Makefile — ising_parallel
#
# Targets
#   make             optimised build (default)
#   make debug       ASan + UBSan build  -> ising_parallel_debug
#   make tsan        ThreadSanitizer     -> ising_parallel_tsan
#   make profile     gprof build         -> ising_parallel_prof
#   make check       headless physics regression tests (check.sh)
#   make install     install binary + man page under $(PREFIX)  (default ~/.local)
#   make uninstall   remove them
#   make clean
#
# Variables you may override on the command line:
#   PREFIX=/usr/local  CC=clang  CFLAGS='-O2'  CPPFLAGS=...  LDFLAGS=...

CC       ?= gcc
PREFIX   ?= $(HOME)/.local
BINDIR    = $(PREFIX)/bin
MANDIR    = $(PREFIX)/share/man/man1

TARGET    = ising_parallel
SRC       = ising_parallel.c
MAN       = ising_parallel.1

# SDL2 flags: sdl2-config if present, else pkg-config
SDL_CFLAGS := $(shell sdl2-config --cflags 2>/dev/null || pkg-config --cflags sdl2 2>/dev/null)
SDL_LIBS   := $(shell sdl2-config --libs   2>/dev/null || pkg-config --libs   sdl2 2>/dev/null)
TTF_CFLAGS := $(shell pkg-config --cflags SDL2_ttf 2>/dev/null)
TTF_LIBS   := $(shell pkg-config --libs   SDL2_ttf 2>/dev/null || echo -lSDL2_ttf)

WARN      = -Wall -Wextra -Wshadow -Wstrict-prototypes -Wno-unused-parameter
CSTD      = -std=c17
CFLAGS   ?= -O3 -march=native -DNDEBUG
ALL_CFLAGS = $(CSTD) $(WARN) $(CFLAGS) $(CPPFLAGS) $(SDL_CFLAGS) $(TTF_CFLAGS)
LIBS      = $(SDL_LIBS) $(TTF_LIBS) -pthread -lm $(LDFLAGS)

.PHONY: all debug tsan profile check install uninstall clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(ALL_CFLAGS) $< -o $@ $(LIBS)

debug: $(SRC)
	$(CC) $(CSTD) $(WARN) -O1 -g3 -fsanitize=address,undefined -fno-omit-frame-pointer \
	    $(CPPFLAGS) $(SDL_CFLAGS) $(TTF_CFLAGS) $< -o $(TARGET)_debug $(LIBS)

tsan: $(SRC)
	$(CC) $(CSTD) $(WARN) -O1 -g3 -fsanitize=thread -fno-omit-frame-pointer \
	    $(CPPFLAGS) $(SDL_CFLAGS) $(TTF_CFLAGS) $< -o $(TARGET)_tsan $(LIBS)

profile: $(SRC)
	$(CC) $(CSTD) $(WARN) -O2 -g -pg $(CPPFLAGS) $(SDL_CFLAGS) $(TTF_CFLAGS) \
	    $< -o $(TARGET)_prof $(LIBS)
	@echo "run ./$(TARGET)_prof --batch ... ; then: gprof $(TARGET)_prof gmon.out | less"

check: $(TARGET)
	./check.sh ./$(TARGET)

install: $(TARGET)
	install -d "$(BINDIR)" "$(MANDIR)"
	install -m 755 $(TARGET) "$(BINDIR)/$(TARGET)"
	install -m 644 $(MAN)    "$(MANDIR)/$(MAN)"
	@echo "installed $(BINDIR)/$(TARGET) and $(MANDIR)/$(MAN)"
	@case ":$$PATH:" in *":$(BINDIR):"*) ;; \
	  *) echo "note: $(BINDIR) is not in PATH";; esac

uninstall:
	rm -f "$(BINDIR)/$(TARGET)" "$(MANDIR)/$(MAN)"
	@echo "removed $(BINDIR)/$(TARGET) and $(MANDIR)/$(MAN)"

clean:
	rm -f $(TARGET) $(TARGET)_debug $(TARGET)_tsan $(TARGET)_prof gmon.out
