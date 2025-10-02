# Makefile for invrate.c on macOS

# Defaults
CC       ?= clang
CSTD     ?= c11
OPTFLAGS ?= -O3 -Ofast -flto -march=native -DNDEBUG
WARN     ?= -Wall -Wextra -Wpedantic

# Runtime params for ./invrate [k] [N] [max_powB] [AFS] [seed]
K        ?= 10
N        ?= 1000000
MAX_POWB ?= 9
AFS      ?= 1
SEED     ?= 12345

# Fallback to Homebrew prefixes
FLINT_PREFIX := $(shell brew --prefix flint 2>/dev/null)
GMP_PREFIX   := $(shell brew --prefix gmp 2>/dev/null)

CFLAGS  += -std=$(CSTD) $(OPTFLAGS) $(WARN) \
			-I$(FLINT_PREFIX)/include -I$(GMP_PREFIX)/include
LDFLAGS += -L$(FLINT_PREFIX)/lib -L$(GMP_PREFIX)/lib
LDLIBS  += -lflint -lgmp

.PHONY: all clean run flags

all: invrate

invrate: invrate.c
	$(CC) $(CFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

run: invrate
	./invrate $(K) $(N) $(MAX_POWB) $(AFS) $(SEED)

flags:
	@echo "CC      = $(CC)"
	@echo "CFLAGS  = $(CFLAGS)"
	@echo "LDFLAGS = $(LDFLAGS)"
	@echo "LDLIBS  = $(LDLIBS)"

clean:
	rm -f invrate *.o *.dSYM