PARA = -std=c++11 -Wall -O3 
com: rwgraph.cpp el2bin.cpp option.cpp GLib.hpp
	g++ -c GLib.hpp -o GLib.o $(PARA)
	g++ -c mappedHeap.hpp -o mappedHeap.o $(PARA)
	g++ -c HeapData.hpp -o HeadData.o $(PARA)
	g++ -c option.cpp -o option.o $(PARA)
	g++ -c rwgraph.cpp -o rwgraph.o $(PARA)
	g++ el2bin.cpp -o el2bin $(PARA)
	g++ verifyInf.cpp rwgraph.o option.o -o verifyInf -fopenmp -std=c++11 -Wall -O3  sfmt/SFMT.c
	g++ verifyDead.cpp rwgraph.o option.o -o verifyDead -fopenmp -std=c++11 -Wall -O3  sfmt/SFMT.c
