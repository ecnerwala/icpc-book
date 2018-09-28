CXX = g++
CXXFLAGS = -Wall -Wextra -pedantic -std=c++11 -O2 -Wshadow -Wformat=2 -Wfloat-equal -Wconversion -Wlogical-op -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align
DEBUGFLAGS = -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fstack-protector -D_FORTIFY_SOURCE=2

CXXFLAGS += $(DEBUGFLAGS)

TARGET := $(notdir $(CURDIR))
EXECUTE := ./$(TARGET)

CASES := $(sort $(basename $(wildcard *.in)))
TESTS := $(sort $(basename $(wildcard *.out)))

all: $(TARGET)

clean:
	-rm -rf $(TARGET) *.res

%: %.cpp
	$(LINK.cpp) $< $(LOADLIBES) $(LDLIBS) -o $@

run: $(TARGET)
	time $(EXECUTE)

%.res: $(TARGET) %.in
	time $(EXECUTE) < $*.in > $*.res

%.out: %

test_%: %.res %.out
	diff $*.res $*.out

runs: $(patsubst %,%.res,$(CASES))
test: $(patsubst %,test_%,$(TESTS))

.PHONY: all clean run test test_% runs

.PRECIOUS: %.res
