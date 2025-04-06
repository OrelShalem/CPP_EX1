# orel8155@gmail.com
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -pedantic
VALGRIND_FLAGS = --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind.log

# Source files
SRCS = graph.cpp UnionFind.cpp Queue.cpp PriorityQueue.cpp
HEADERS = graph.hpp UnionFind.hpp Queue.hpp PriorityQueue.hpp

# Targets
MAIN_TARGET = main
TEST_TARGET = test

# Default target
all: $(MAIN_TARGET)

# Compile and run the demonstration file
main: $(MAIN_TARGET).cpp $(SRCS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(MAIN_TARGET) $(MAIN_TARGET).cpp $(SRCS)
	./$(MAIN_TARGET)

# Compile and run the test file
test: $(TEST_TARGET).cpp $(SRCS) $(HEADERS) doctest.h
	$(CXX) $(CXXFLAGS) -o $(TEST_TARGET) $(TEST_TARGET).cpp $(SRCS)
	./$(TEST_TARGET)

# Check for memory leaks with valgrind
valgrind: $(MAIN_TARGET).cpp $(SRCS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(MAIN_TARGET) $(MAIN_TARGET).cpp $(SRCS)
	valgrind $(VALGRIND_FLAGS) ./$(MAIN_TARGET)

# Check for memory leaks in tests
valgrind-test: $(TEST_TARGET)
	valgrind $(VALGRIND_FLAGS) ./$(TEST_TARGET)

# Clean the project
clean:
	rm -f $(MAIN_TARGET) $(TEST_TARGET) *.o

.PHONY: all main test valgrind valgrind-test clean
