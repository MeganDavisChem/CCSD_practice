#include "hf.h"
#include "molecule.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

int main()
{
	//molecule name, basis set, nelec, e conv, rms conv
	hf fock("h2o", "STO-3G", 10, 1e-10, 1e-10);
	fock.fock_procedure(); 
//	fock.mp2_proc();
	//Econv and RMS conv
//	fock.cc_proc(1e-10,1e-10);

	//(T)
//	fock.full_pt();

	//MO transformation
//	fock.mo_fock();
//	fock.smarty();
	fock.noddy();
	fock.spin_orbit_CI();
	fock.spin_mo_fock();

//	fock.CI();
	fock.CI_spin();


	return 0;
}
