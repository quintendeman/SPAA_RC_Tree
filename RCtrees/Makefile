# Variables
CXX = g++
CXXFLAGS = -O3 -std=c++17 -pthread -g
INCLUDES = -I ../include/
TARGET = RC
SRC = RC.cpp

# Default target
all: $(TARGET)

# Rule to build the target
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET) $(SRC)

# Clean up build files
clean:
	rm -f $(TARGET)

# Phony targets
.PHONY: all clean
