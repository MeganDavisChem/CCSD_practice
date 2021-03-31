#include "hf.h"
#include "molecule.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

int main()
{
	//molecule name, basis set, nelec, e conv, rms conv
	hf fock("h2o", "DZP", 10, 1e-10, 1e-10);
	fock.fock_procedure();
	fock.mp2_proc();

	//Econv and RMS conv
	fock.cc_proc(1e-10,1e-10);

	//(T)
	fock.full_pt();

	return 0;
}
