harmonic : harmonic.o molecule.o 
	  c++ -o harmonic molecule.o harmonic.o

molecule.o : molecule.cc molecule.h masses.h
			c++ -c molecule.cc
harmonic.o : harmonic.cc molecule.h masses.h
	   c++ -c harmonic.cc
