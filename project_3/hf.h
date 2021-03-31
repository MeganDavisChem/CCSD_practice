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

		//Creation and annihiliation
		hf(string mol, string basis, int elec, double e_conv_crit_input, double rms_conv_crit_input);
		~hf();
};
