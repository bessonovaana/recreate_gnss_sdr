CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -Imath -Iacquistion -Itracking

MATH_SRC = math/correlation.cpp math/generate_code.cpp math/iq_samples.cpp
COMMON_SRC = acquistion/1.cpp tracking/1.cpp $(MATH_SRC)

.PHONY: all clean

all: main

main: main.cpp $(COMMON_SRC)
	$(CXX) $(CXXFLAGS) main.cpp $(COMMON_SRC) -o main

clean:
	rm -f main acquistion/acquisition tracking/tracking
