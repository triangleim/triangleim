PARA = -std=c++11 -Wall -O3
com: rwgraph.cpp SSA.cpp el2bin.cpp option.cpp GLib.hpp Sandwich.cpp
	g++ -c GLib.hpp -o GLib.o $(PARA)
	g++ -c mappedHeap.hpp -o mappedHeap.o $(PARA)
	g++ -c HeapData.hpp -o HeadData.o $(PARA)
	g++ -c option.cpp -o option.o $(PARA)
	g++ -c rwgraph.cpp -o rwgraph.o $(PARA)
	# # g++ SSA.cpp rwgraph.o option.o -o SSA -fopenmp $(PARA) sfmt/SFMT.c
	# g++ verifyInf.cpp rwgraph.o option.o -o verifyInf -fopenmp $(PARA) sfmt/SFMT.c
	# g++ verifyDead.cpp rwgraph.o option.o -o verifyDead -fopenmp $(PARA) sfmt/SFMT.c
	# g++ el2bin.cpp -o el2bin $(PARA)
	# g++ format.cpp -o format $(PARA)
	g++ Sandwich.cpp rwgraph.o option.o -o Sandwich -fopenmp $(PARA) sfmt/SFMT.c

