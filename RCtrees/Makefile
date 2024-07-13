# Variables
CXX = g++
CXXFLAGS = -O3 -std=c++17 -pthread -g
INCLUDES = -I ../include/ -I /usr/local/include/jemalloc/
LDFLAGS = -ljemalloc
TARGET = RC.out
SRC = RC.cpp

# Default target
all: $(TARGET)

# Rule to build the target
$(TARGET): $(SRC) RC.h cluster.h RCdynamic.h adjacency_linked_list.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET) $(SRC)

# Clean up build files
clean:
	rm -f $(TARGET)

# Phony targets
.PHONY: all clean
