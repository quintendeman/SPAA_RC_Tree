# Variables
CXX = g++
CXXFLAGS = -O3 -std=c++17 -pthread -g -Werror -Wparentheses -flto
# CXXFLAGS = -O0 -flto -std=c++17 -pthread -g -Werror -Wparentheses
INCLUDES = -I ../include/ -I /usr/local/include/jemalloc/ -I parhash_include/ -I LCA -I .
LDFLAGS = -ljemalloc
TARGET = RC.out
SRC = RC.cpp

# Default target
# all: RC.out testMST.out  testHashing.out testTreeGen.out testBatchedSubtree.out testPathQueries.out testBatchInsertion.out testMST.out testMSTfromFile.out testNewTernarizer.out
all: testMST.out   testTreeGen.out testBatchedSubtree.out testPathQueries.out testBatchInsertion.out testMST.out lca.out testMSTfromFile.out testMSTfromFile2.out

# Rule to build the target
# RC.out: RC.cpp RC.h cluster.h RCdynamic.h adjacency_linked_list.h incMST.h random_trees.h utils.h
# 	$(CXX) $(CXXFLAGS) $(INCLUDES) -o RC.out RC.cpp

testMST.out: testMST.cpp ternarizer.h treeGen.h incMST.h RC.h cluster.h RCdynamic.h  adjacency_linked_list.h  random_trees.h utils.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o testMST.out testMST.cpp

testMSTfromFile.out: testMSTfromFile.cpp ternarizer.h treeGen.h incMST.h RC.h cluster.h RCdynamic.h  adjacency_linked_list.h  random_trees.h utils.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o testMSTfromFile.out testMSTfromFile.cpp

testMSTfromFile2.out: testMSTfromFile2.cpp ternarizer.h treeGen.h incMST.h RC.h cluster.h RCdynamic.h  adjacency_linked_list.h  random_trees.h utils.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o testMSTfromFile2.out testMSTfromFile2.cpp


testTreeGen.out: testTreeGen.cpp treeGen.h RC.h cluster.h RCdynamic.h  adjacency_linked_list.h  random_trees.h utils.h
	$(CXX)  $(CXXFLAGS) $(INCLUDES) -o testTreeGen.out testTreeGen.cpp

# testHashing.out: treeGen.h RC.h cluster.h RCdynamic.h adjacency_linked_list.h  random_trees.h utils.h testHashing.cpp ternarizer.h
# 	$(CXX)  $(CXXFLAGS)  -o testHashing.out testHashing.cpp $(INCLUDES)

testBatchedSubtree.out: treeGen.h RC.h cluster.h RCdynamic.h adjacency_linked_list.h  random_trees.h utils.h testBatchedSubtree.cpp ternarizer.h subtree_query.h
	$(CXX)  $(CXXFLAGS)  -o testBatchedSubtree.out testBatchedSubtree.cpp $(INCLUDES)

testPathQueries.out: testPathQueries.cpp treeGen.h RC.h cluster.h RCdynamic.h adjacency_linked_list.h  random_trees.h utils.h ternarizer.h 
	$(CXX)  $(CXXFLAGS)  -o testPathQueries.out testPathQueries.cpp $(INCLUDES)

testBatchInsertion.out: treeGen.h RC.h cluster.h RCdynamic.h adjacency_linked_list.h  random_trees.h utils.h subtree_query.h testBatchInsertion.cpp
	$(CXX)  $(CXXFLAGS)  -o testBatchInsertion.out testBatchInsertion.cpp $(INCLUDES)


testNewTernarizer.out: testNewTernarizer.cpp treeGen.h RC.h cluster.h RCdynamic.h adjacency_linked_list.h  random_trees.h utils.h ternarizer.h incMST.h
	$(CXX)  $(CXXFLAGS)  -o testNewTernarizer.out testNewTernarizer.cpp $(INCLUDES)

RC_CREATE_DEPS = MIS.h RC.h RCdynamic.h RC_test.h RCdynamic.h
RC_BASE_DEPS = cluster.h adjacency_linked_list.h 
UTIL_DEPS = random_trees.h utils.h treeGen.h treePrimitives.h
BASE_DEPS = $(RC_CREATE_DEPS) $(RC_BASE_DEPS) $(UTIL_DEPS)
QUERY_DEPS = subtree_query.h path_query.h 
MST_DEPS = incMST.h
LCA_DEPS = LCA/VanillaLCA.h LCA/static_seqn_LCA.h LCA/static_par_LCA.h LCA/LCA.h LCA/fixedLCA.h LCA/LCA_test.h

# Default target

splca.out : static_par_LCA.cpp $(LCA_DEPS) $(BASE_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o splca.out static_par_LCA.cpp

sslca.out : static_seqn_LCA.cpp $(LCA_DEPS) $(BASE_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o sslca.out static_seqn_LCA.cpp

#Linux case sensitive, LCA.cpp not lca.cpp
lca.out : LCA.cpp $(LCA_DEPS) $(BASE_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o lca.out LCA.cpp

# Clean up build files
clean:
	rm -f *.out

# Phony targets
.PHONY: all clean
