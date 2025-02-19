# Variables
CXX = g++
CXXFLAGS = -O3 -std=c++17 -pthread -g -Werror -Wparentheses
# CXXFLAGS = -O0 -flto -std=c++17 -pthread -g -Werror -Wparentheses
INCLUDES = -I ../include/ -I /usr/local/include/jemalloc/
LDFLAGS = -ljemalloc
TARGET = RC.out
SRC = RC.cpp

# Default target
all: RC.out testMST.out  testHashing.out testTreeGen.out

# Rule to build the target
RC.out: RC.cpp RC.h cluster.h RCdynamic.h adjacency_linked_list.h incMST.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o RC.out RC.cpp

testMST.out: testMST.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o testMST.out testMST.cpp

testTreeGen.out: testTreeGen.cpp treeGen.h RC.h cluster.h RCdynamic.h adjacency_linked_list.h 
	$(CXX)  $(CXXFLAGS) $(INCLUDES) -o testTreeGen.out testTreeGen.cpp

testHashing.out: treeGen.h RC.h cluster.h RCdynamic.h adjacency_linked_list.h testHashing.cpp ternarizer.h
	$(CXX)  $(CXXFLAGS)  -o testHashing.out testHashing.cpp $(INCLUDES)


# Clean up build files
clean:
	rm -f *.out

# Phony targets
.PHONY: all clean
