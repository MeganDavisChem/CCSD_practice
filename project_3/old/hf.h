#include <string>
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Eigenvalues"
#include "eigen3/Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;


using namespace std;

class hf
{
	public:
		//Variables
		double enuc;
		Matrix ess;
		double **S;
		double **T;
		double **V;
		double **Hcore;

		//Two e integrals very unsure about this
		double ****eri;
		
		//Function declarations
		void example();
		
		//Creation and annhiliation
		hf(const char *mol, const char *basis, int q);
		~hf();
};

