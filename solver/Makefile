# Compiler and flags
FC = mpif90
# A100
# FLAGS = -fast -lm -acc -gpu=cc80,cuda11.2,fastmath -Minfo=accel
# GH200
FLAGS = -fast -lm -acc -gpu=cc90,cuda12.3,fastmath -Minfo=accel

# Source and object files
SRC = $(wildcard src/*.f90)
OBJ = $(SRC:src/%.f90=build/%.o)
# OBJ = build/solver.o build/main.o

# Output executable
TARGET = main

# Default target
all: $(TARGET)

# Rule to build the executable
$(TARGET): $(OBJ)
	$(FC) $(FLAGS) $^ -o $@

# Rule to create object files with explicit dependency
build/solver.o: src/solver.f90 | build
	$(FC) $(FLAGS) -c $< -o $@

build/solver_acc.o: src/solver_acc.f90 | build
	$(FC) $(FLAGS) -c $< -o $@

build/main.o: src/main.f90 build/solver.o build/solver_acc.o | build
	$(FC) $(FLAGS) -c $< -o $@

# Rule to create the build directory if it doesn't exist
build:
	mkdir -p build

# Clean target to remove the generated files
clean:
	rm -rf build $(TARGET) *.mod src/*.mod

.PHONY: all clean build