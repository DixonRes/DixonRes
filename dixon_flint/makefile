# Compiler
CC = gcc

# Compiler flags
CFLAGS = -O3 -march=native -fopenmp

# Linker flags
LDFLAGS = -fopenmp

# Libraries
LIBS = -lflint -lmpfr -lgmp -lpthread -lstdc++ -lpml2

# Custom library path
LIBPATH = -L/home/suohaohai02/mylinks

# Source files
DIXON_SRC = dixon.c
RESCUE_SRC = rescue_attack_dixon.c

# Output executables
DIXON_TARGET = dixon
RESCUE_TARGET = rescue

# Default target (run with just 'make')
default: $(DIXON_TARGET)

# Build both targets
all: $(DIXON_TARGET) $(RESCUE_TARGET)

$(DIXON_TARGET): $(DIXON_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBPATH) $(LIBS) $(LDFLAGS)

$(RESCUE_TARGET): $(RESCUE_SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBPATH) $(LIBS) $(LDFLAGS)

clean:
	rm -f $(DIXON_TARGET) $(RESCUE_TARGET)

.PHONY: default all clean