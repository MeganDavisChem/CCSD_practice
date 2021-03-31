#include "hf.h"
#include "molecule.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

int main()
{
	//molecule name, basis set, nelec, e conv, rms conv
	hf fock("h2o", "STO-3G", 10, 1e-12, 1e-10);
	fock.fock_procedure();
	fock.mp2_proc();
	return 0;
}
