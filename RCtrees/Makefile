# Variables
CXX = g++
# CXXFLAGS = -O0 -flto -std=c++17 -pthread -g -Werror -Wparentheses
CXXFLAGS = -O0 -std=c++17 -pthread -g -Werror -Wparentheses

#must include -I (includes) before each new item!
INCLUDES = -I ../include/ -I /usr/local/include/jemalloc/ -I include/ 
LDFLAGS = -ljemalloc
TARGET = RC.out
SRC = RC.cpp

TARGET2 = LCA.out
SRC2 = LCA.cpp
RC_DEPS = RC.h cluster.h RCdynamic.h adjacency_linked_list.h RC_test.h 
LCA_DEPS = LCA.h VanillaLCA.h static_seqn_LCA.h static_par_LCA.h random_trees.h

# Default target
all: RC.out testMST.out testTreeGen.out

lca: $(TARGET2)

splca : splca.out

sslca : sslca.out

splca.out : static_par_LCA.cpp static_par_LCA.h static_seqn_LCA.h VanillaLCA.h RC.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o splca.out static_par_LCA.cpp


sslca.out : static_seqn_LCA.cpp static_seqn_LCA.h RC.h VanillaLCA.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o sslca.out static_seqn_LCA.cpp

$(TARGET2): $(SRC2)  $(RC_DEPS) $(LCA_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET2) $(SRC2)

# Rule to build the target
RC.out: RC.cpp RC.h cluster.h RCdynamic.h adjacency_linked_list.h incMST.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o RC.out RC.cpp

testMST.out: testMST.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o testMST.out testMST.cpp

testTreeGen.out: testTreeGen.cpp treeGen.h RC.h cluster.h RCdynamic.h adjacency_linked_list.h 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o testTreeGen.out testTreeGen.cpp


# Clean up build files
clean:
	rm -f RC.out testMST.out testTreeGen.out
	rm -f $(TARGET) 
	rm -f $(TARGET2)
	rm -f sslca.out

# Phony targets
.PHONY: all clean
