# Compiler
CC = gcc

# Compiler flags  -fsanitize=address -static-libasan
CFLAGS = -O3 -march=native -fopenmp

# Linker flags
LDFLAGS = -fopenmp

# Libraries
LIBS = -lflint -lmpfr -lgmp -lpthread -lstdc++ -lpml2

# Custom library path
LIBPATH = -L/home/suohaohai02/mylinks

# Source files
SRC = dixon.c

# Output executable
TARGET = dixon

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBPATH) $(LIBS) $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: all clean
