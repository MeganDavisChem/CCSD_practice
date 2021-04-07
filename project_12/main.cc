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
	fock.smarty();
//	fock.noddy();

//	Pick spin orbits
//	fock.spin_orbit();



	//Pick CI (match w/ orbits)
//	fock.CI();

//	fock.spin_orbit();
//	fock.spin_mo_fock();

	fock.spin_orbit_CI();
//	fock.spin_orbit_CI();
	fock.CI_spin();

	//mo_fock after spin orbit
//


//fock.spin_orbit();

	//RPA
//	fock.RPA();
	fock.RPA_reduced();


	return 0;
}
