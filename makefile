CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic
VALGRIND_FLAGS = --leak-check=full --show-leak-kinds=all --track-origins=yes

# קבצי מקור
SRCS = graph.cpp
HEADERS = graph.hpp

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
test: $(TEST_TARGET).cpp $(SRCS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(TEST_TARGET) $(TEST_TARGET).cpp $(SRCS)
	./$(TEST_TARGET)

# בדיקת זליגת זיכרון עם valgrind
valgrind: $(MAIN_TARGET)
	valgrind $(VALGRIND_FLAGS) ./$(MAIN_TARGET)

# ניקוי הפרויקט
clean:
	rm -f $(MAIN_TARGET) $(TEST_TARGET) *.o

.PHONY: all main test valgrind clean
