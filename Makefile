all: gc
gc: Main.cpp
	g++ -std=c++11 Main.cpp -o gcp.out
clean:
	rm -f gc