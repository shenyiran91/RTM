# Files
EXEC := simulation
SRC	 := $(wildcard *.cpp)
OBJ  := $(patsubst %.cpp, %.o, $(SRC))

# Options
CXX      := g++
CPPFLAGS := -Wall -g -O3 -std=c++17
INCLUDES := -I/usr/include/eigen3
LDLIBS   := -lm

# Rules
$(EXEC): $(OBJ)
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(LDLIBS) -o $@ $^
%.o: %.c
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(LDLIBS) -c $<

.PHONY: clean
clean:
	$(RM) $(OBJ) $(EXEC)
