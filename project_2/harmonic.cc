#include "molecule.h"
#include <iostream>
#include <cmath>
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Eigenvalues"
#include "eigen3/Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

using namespace std;

int main()
{
	std::cout.precision(7);
	std::cout.setf(std::ios::fixed);
	Molecule mol("geom.dat",0);
	cout << "geo\n" << endl;
	mol.print_geom();

	Matrix Hess(mol.natom*3,mol.natom*3);
	printf("\n");
	mol.mass_weight();
	for (int i=0; i < mol.natom*3; i++)
	{
		for (int j=0; j < mol.natom*3; j++)
		{
			Hess(i,j) = mol.H[i][j];
		}
	}
	Eigen::SelfAdjointEigenSolver<Matrix> solver(Hess);
	Matrix evecs = solver.eigenvectors();
	Matrix evals = solver.eigenvalues();
	cout << "mass weighted Hessian \n" << endl;
	cout << Hess << endl;
	cout << "Hessian eigenvalues \n" << endl;
	cout << evals << endl;
	cout << "Frequencies in cm-1 \n" << endl;
	cout << evals.array().sqrt()*5140.48400265 << endl;
	return 0;
}
