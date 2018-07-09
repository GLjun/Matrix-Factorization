CXX = g++-7
CXXFLAGS = -O3 -std=c++11 -g

all : main

main.o : matrix_factorization.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $^

main : main.o
	$(CXX) $(CXXFLAGS) -o $@ $<

clean :
	rm -f *.o ./main

