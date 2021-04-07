#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Eigenvalues"
#include "eigen3/Eigen/Core"
#include <string>
#define INDEX(i,j) (i>j) ? (ioff[i]+j) : (ioff[j]+i)


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;


using namespace std;

class hf
{
	public:
		int orbitals;
		int nelec;

		double e_conv_crit;
		double rms_conv_crit;

		double enuc;
		double **S;
		double **T;
		double **V;
		double *TEI;

		double **dip_x;
		double **dip_y;
		double **dip_z;


		double e_elec;
		double e_total;
		double e_elec_prev;
		double e_total_prev;
		double deltaE;
		double rms;

		Matrix S_12;
		Matrix F_prime;
		Matrix Hcore;
		Matrix Density;
		Matrix F;
		Matrix D_old;
		Matrix C;
		Matrix mo_F;
		int *ioff;

		//MP2 variables
		double *TEI_MO;
		double Emp2;
		
		//CCSD variables
		double ****ints;
		int nmo;
		double **smo_F;
		double **T1;
		double **T1_old;
		double ****T2;
		double ****T2_old;
		
		//Stanton paper intermediates
		double **Fae;
		double **Fmi;
		double **Fme;
		double ****Wmnij;
		double ****Wabef;
		double ****Wmbej;

		//Stanton tau
		double ****tilde_tau;
		double ****tau;
		double **Dia;
		double ****Dijab;

		double cc_e;
		double cc_e_old;
		int iter;

		//(T) stuff
//		double ******Dijkabc;
//		double ******T3D;
//		double ******T3C;
		double e_pt;

		//CI variables
		Matrix HCI;
		Matrix HCI_spin;
		Matrix HRPA;
		Matrix HRPA_reduced;
	

		//Individual functions
		void ortho();	
		void guess_density();
		void scf_energy();
		void new_fock();
		void mo_fock();
		void fock_procedure();
		double dipole(int cart);
		bool e_conv_check(double e_conv);
		bool rms_conv_check(double rms_conv);


		//MP2 part (lazy)
		void noddy();
		void smarty();
		void mp2_energy();
		void mp2_proc();

		//CCSD
		void spin_orbit();
		void spin_orbit_CI();
		void spin_mo_fock();

		void init_amp();
		void cc_inter();
		void amp_up();
		void cc_energy();
		void test_mo_fock();
		void cc_iterate();
		bool cc_conv_check(double cc_rms_conv);
		void cc_proc(double ccEconv,double ccRMSconv);

		//(T)
		void full_pt();

		//CI
		void CI();
		void CI_spin();

		//RPA
		void RPA();
		void RPA_reduced();


		//Creation and annihiliation
		hf(string mol, string basis, int elec, double e_conv_crit_input, double rms_conv_crit_input);
		~hf();
};
