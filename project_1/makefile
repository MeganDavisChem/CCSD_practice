geo : geo.o molecule.o 
	  c++ -o geo molecule.o geo.o

molecule.o : molecule.cc molecule.h masses.h
			c++ -c molecule.cc
geo.o : geo.cc molecule.h masses.h
	   c++ -c geo.cc
