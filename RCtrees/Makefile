# Variables
CXX = g++
# CXXFLAGS = -O0 -flto -std=c++17 -pthread -g -Werror -Wparentheses
CXXFLAGS = -O0 -std=c++17 -pthread -g -Werror -Wparentheses

#must include -I (includes) before each new item!
INCLUDES = -I ../include/ -I /usr/local/include/jemalloc/ -I parhash_include/ -I LCA -I .
LDFLAGS = -ljemalloc

RC_CREATE_DEPS = MIS.h RC.h RCdynamic.h RC_test.h RCdynamic.h
RC_BASE_DEPS = cluster.h adjacency_linked_list.h 
UTIL_DEPS = random_trees.h utils.h treeGen.h
BASE_DEPS = $(RC_CREATE_DEPS) $(RC_BASE_DEPS) $(UTIL_DEPS)
QUERY_DEPS = subtree_query.h path_query.h 
MST_DEPS = incMST.h
LCA_DEPS = LCA/VanillaLCA.h LCA/static_seqn_LCA.h LCA/static_par_LCA.h LCA/LCA.h LCA/fixedLCA.h

# Default target
all: RC.out testMST.out testTreeGen.out splca.out sslca.out lca.out

splca.out : static_par_LCA.cpp $(LCA_DEPS) $(BASE_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o splca.out static_par_LCA.cpp

sslca.out : static_seqn_LCA.cpp $(LCA_DEPS) $(BASE_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o sslca.out static_seqn_LCA.cpp

#Linux case sensitive, LCA.cpp not lca.cpp
lca.out : LCA.cpp $(LCA_DEPS) $(BASE_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o lca.out LCA.cpp


RC.out: RC.cpp $(BASE_DEPS) $(MST_DEPS) $(QUERY_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o RC.out RC.cpp

testMST.out: testMST.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o testMST.out testMST.cpp

testTreeGen.out: testTreeGen.cpp $(BASE_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o testTreeGen.out testTreeGen.cpp


# Clean up build files
clean:
	rm -f RC.out testMST.out testTreeGen.out sslca.out splca.out lca.out
	

# Phony targets
.PHONY: all clean
