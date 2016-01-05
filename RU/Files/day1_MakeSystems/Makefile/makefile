CXX  = g++ -c -I. -O3 -m64
LINK = g++ -lm 

main: square.obj circle.obj main.obj
	$(LINK) main.obj square.obj circle.obj -omain

circle.obj: circle.cxx fprotos.hxx
	$(CXX) circle.cxx -ocircle.obj

square.obj: square.cxx fprotos.hxx
	$(CXX) square.cxx -osquare.obj

main.obj: main.cxx fprotos.hxx
	$(CXX) main.cxx -omain.obj


all: main


clean:
	rm -rf *.obj
	rm -rf main

