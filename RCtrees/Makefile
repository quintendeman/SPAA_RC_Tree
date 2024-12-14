# Variables
CXX = g++
CXXFLAGS = -O0 -std=c++17 -pthread -g -Werror -Wparentheses
INCLUDES = -I ../include/ -I /usr/local/include/jemalloc/
LDFLAGS = -ljemalloc
TARGET = RC.out
SRC = RC.cpp

TARGET2 = LCA.out
SRC2 = LCA.cpp
RC_DEPS = RC.h cluster.h RCdynamic.h adjacency_linked_list.h 
LCA_DEPS = LCA.h VanillaLCA.h 

# Default target
all: $(TARGET)

lca: $(TARGET2)

sslca : sslca.out

sslca.out : static_seqn_LCA.cpp static_seqn_LCA.h RC.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o sslca.out static_seqn_LCA.cpp

$(TARGET2): $(SRC2)  $(RC_DEPS) $(LCA_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET2) $(SRC2)

# Rule to build the target
$(TARGET): $(SRC) RC.h cluster.h RCdynamic.h adjacency_linked_list.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET) $(SRC)

# Clean up build files
clean:
	rm -f $(TARGET) 
	rm -f $(TARGET2)
	rm -f sslca.out

# Phony targets
.PHONY: all clean
