#Created by Rui 11/29/2020

CC = c++
CPPFLAGS =-g -Wall -std=c++14 -O3
INCLUDES = -I/usr/local/include -L/usr/local/lib -lboost_unit_test_framework -static -lpthread
LIB = -I/usr/local/include -L/usr/local/lib -lpthread

TEST_CASES := binpacking_test search_test io_test glycan_builder_test glycan_test protein_test  
TEST_CASES_2 := precursor_match_test search_sequence_test search_glycan_test search_engine_test multi_comparison_test
TEST_CASES_3 := lsh_clustering_test


search:
	$(CC) $(CPPFLAGS) -o searching \
	apps/searching.cpp model/glycan/nglycan_complex.cpp model/glycan/nglycan_hybrid.cpp $(LIB)

search_windows:
	$(CC) $(CPPFLAGS) -o searching_windows \
	apps/searching_windows.cpp model/glycan/nglycan_complex.cpp model/glycan/nglycan_hybrid.cpp $(LIB)

binpacking_test:
	$(CC) $(CPPFLAGS) -o test/binpacking_test \
	 engine/spectrum/binpacking_test.cpp $(INCLUDES)

search_test:
	$(CC) $(CPPFLAGS) -o test/search_test \
	algorithm/search/search_test.cpp $(INCLUDES)

io_test:
	$(CC) $(CPPFLAGS) -o test/io_test \
	util/io/io_test.cpp  $(INCLUDES)

glycan_builder_test:
	$(CC) $(CPPFLAGS) -o test/glycan_builder_test \
	engine/glycan/builder_test.cpp model/glycan/nglycan_complex.cpp $(INCLUDES)

glycan_test:
	$(CC) $(CPPFLAGS) -o test/glycan_test \
	 model/glycan/glycan_test.cpp model/glycan/nglycan_complex.cpp model/glycan/nglycan_hybrid.cpp $(INCLUDES)

protein_test:
	$(CC) $(CPPFLAGS) -o test/protein_test \
	engine/protein/protein_test.cpp  $(INCLUDES)

precursor_match_test:
	$(CC) $(CPPFLAGS) -o test/precursor_match_test \
	engine/search/precursor_match_test.cpp model/glycan/nglycan_complex.cpp $(INCLUDES)

search_sequence_test:
	$(CC) $(CPPFLAGS) -o test/search_sequence_test \
	engine/search/search_sequence_test.cpp model/glycan/nglycan_complex.cpp $(INCLUDES)

search_glycan_test:
	$(CC) $(CPPFLAGS) -o test/search_glycan_test \
	engine/search/search_glycan_test.cpp model/glycan/nglycan_complex.cpp $(INCLUDES)

search_engine_test:
	$(CC) $(CPPFLAGS) -o test/search_engine_test \
	engine/search/search_engine_test.cpp model/glycan/nglycan_complex.cpp model/glycan/nglycan_hybrid.cpp $(INCLUDES)

multi_comparison_test:
	$(CC) $(CPPFLAGS) -o test/multi_comparison_test \
	engine/analysis/multi_comparison_test.cpp $(INCLUDES)

lsh_clustering_test:
	$(CC) $(CPPFLAGS) -o test/lsh_clustering_test \
	engine/spectrum/lsh_clustering_test.cpp $(INCLUDES)

modification_test:
	$(CC) $(CPPFLAGS) -o test/modification_test \
	engine/protein/modification_test.cpp $(INCLUDES)

# test
test: ${TEST_CASES} ${TEST_CASES_2} ${TEST_CASES_3}

# clean up
clean:
	rm -f core test/* *.o clustering searching