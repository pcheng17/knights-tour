# Makefile for Knight's Tour

# Compiler selection - tries clang first, falls back to gcc
CC := $(shell which clang 2>/dev/null || which gcc)

# Compiler flags
CFLAGS = -Wall -Wextra -std=c99 -O3

# Source and target
SRC = knights_tour.c
TARGET = knights_tour

# Default target
all: $(TARGET)

# Build rule
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

# Clean rule
clean:
	rm -f $(TARGET)

# Force rebuild
rebuild: clean all

# Phony targets
.PHONY: all clean rebuild
