#Created by Rui 5/17/20

CC = c++
CPPFLAGS =-g -Wall -std=c++14 -O3
INCLUDES = -I/usr/local/include -L/usr/local/lib -lboost_unit_test_framework -static -lpthread
LIB = -I/usr/local/include -L/usr/local/lib -lpthread

TEST_CASES := binpacking_test search_test io_test

binpacking_test:
	$(CC) $(CPPFLAGS) -o test/binpacking_test \
	 engine/spectrum/binpacking_test.cpp $(INCLUDES)

search_test:
	$(CC) $(CPPFLAGS) -o test/search_test \
	algorithm/search/search_test.cpp $(INCLUDES)

io_test:
	$(CC) $(CPPFLAGS) -o test/io_test \
	util/io/io_test.cpp  $(INCLUDES)

# test
test: ${TEST_CASES}

# clean up
clean:
	rm -f core test/* *.o clustering searching