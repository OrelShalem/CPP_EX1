CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic
VALGRIND_FLAGS = --leak-check=full --show-leak-kinds=all --track-origins=yes

# קבצי מקור
SRCS = graph.cpp UnionFind.cpp Queue.cpp PriorityQueue.cpp
HEADERS = graph.hpp UnionFind.hpp Queue.hpp PriorityQueue.hpp

# יעדים
MAIN_TARGET = main
TEST_TARGET = test

# פקודת ברירת מחדל
all: $(MAIN_TARGET)

# קומפילציה והרצה של קובץ ההדגמה
main: $(MAIN_TARGET).cpp $(SRCS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(MAIN_TARGET) $(MAIN_TARGET).cpp $(SRCS)
	./$(MAIN_TARGET)

# קומפילציה והרצה של קובץ הבדיקות
test: $(TEST_TARGET).cpp $(SRCS) $(HEADERS) doctest.h
	$(CXX) $(CXXFLAGS) -o $(TEST_TARGET) $(TEST_TARGET).cpp $(SRCS)
	./$(TEST_TARGET)

# בדיקת זליגת זיכרון עם valgrind
valgrind: $(MAIN_TARGET)
	valgrind $(VALGRIND_FLAGS) ./$(MAIN_TARGET)

# בדיקת זליגת זיכרון בטסטים
valgrind-test: $(TEST_TARGET)
	valgrind $(VALGRIND_FLAGS) ./$(TEST_TARGET)

# ניקוי הפרויקט
clean:
	rm -f $(MAIN_TARGET) $(TEST_TARGET) *.o

.PHONY: all main test valgrind valgrind-test clean
