#include "hf.h"
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
	hf fock("h2o","STO-3G",0);

	fock.ess.resize(3,3);

	for(int i=0; i < 3; i++)
	{
		for(int j=0; j < 3; j++)
			fock.ess(i,j) = i + j;
	}

	cout << fock.ess << endl;
}
