hf : main.o hf.o molecule.o
	c++ -o hf hf.o main.o molecule.o

main.o : main.cc hf.h masses.h
	c++ -c main.cc

hf.o : hf.cc hf.h masses.h
	c++ -c hf.cc

molecule.o : molecule.cc molecule.h masses.h
	c++ -c molecule.cc

