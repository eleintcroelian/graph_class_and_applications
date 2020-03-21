# conventional variable for c++ compiler
CXX := clang++

# conventional variable for C preprocessor
CPPFLAGS := 

# conventional variable for C++ compiler flags
CXXFLAGS := -O3 -std=c++11 -Wall -Wextra -Wconversion -pedantic

# conventional variable for linker flags
LDFLAGS := -ljpeg

TARGET := main
OBJS := main.o
INCS := MinStack.hpp

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LDFLAGS)

%.o: %.cpp $(INCS)
	$(CXX) -c -o $@ $< $(CPPFLAGS) $(CXXFLAGS)

# use .PHONY for targets that do not produce a file
.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET) *~

.PHONY: run
run: $(TARGET)
	./$(TARGET)