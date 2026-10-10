# Compiler and flags
CC      = gcc
CFLAGS  = -std=c17 -Wall -Wextra -pedantic -g
LDLIBS  = -lm

# Project files
TARGET  = mySolver
SRCS    = mySolver.c solverType.c
HDRS    = solverType.h
OBJS    = $(SRCS:.c=.o)

# These targets are commands, not files
.PHONY: all clean

# Default target: runs when you just type "make"
all: $(TARGET)

# Link all object files into the final program
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

# Compile each .c into a .o; rebuild if any header changes
%.o: %.c $(HDRS)
	$(CC) $(CFLAGS) -c $< -o $@

# Delete build output
clean:
	rm -f $(TARGET) $(OBJS)
