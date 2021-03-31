fock : fock.o hf.o 
	  c++ -o fock hf.o fock.o

hf.o : hf.cc hf.h masses.h
			c++ -c hf.cc
fock.o : fock.cc hf.h masses.h
	   c++ -c fock.cc
