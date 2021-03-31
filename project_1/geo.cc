#include "molecule.h"
#include <iostream>
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Eigenvalues"
#include "eigen3/Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

using namespace std;

void hello()
{
	cout << "hello world!" << endl;

}

int main(int argc, char *argv[])
{
	Molecule mol("geom.dat",0);
	cout << "Number of atoms: " << mol.natom << endl;
	cout << "Input cartesian geometry:" << endl;
	mol.print_geom();	

	cout << "Bond distances: " << endl;
	for(int i=0; i < mol.natom; i++)
	{
		for(int j=0; j < i; j++)
		{
			mol.bonds[i][j] = mol.bond_lengths(i,j);
			cout << i << "," << j << ": " << mol.bonds[i][j] << endl;
		}
	}
	printf("bond angles: \n");
	for(int i=0; i < mol.natom; i++) 
	{
		for(int j=0; j < i; j++)
		{
			for(int k=0; k < j; k++)
			{
				if(mol.bonds[i][j] < 4.0 && mol.bonds[j][k] < 4.0)
					printf("%2d-%2d-%2d %10.6f\n", i, j, k, mol.bond_angles(i,j,k));
			}
		}
	}

	printf("out-of-plane angles: \n");
	for(int i=0; i < mol.natom; i++)
	{
		for(int k=0; k < mol.natom; k++)
		{
			for(int j=0; j < mol.natom; j++)
			{
				for(int l=0; l < j; l++)
					if(i != j && i !=k && i!=l && j!=k && k!=l && mol.bond_lengths(i,k) < 4.0 && mol.bond_lengths(k,j) < 4.0 && mol.bond_lengths(k,l) < 4.0)
						printf("%2d-%2d-%2d-%2d %10.6f\n",i,j,k,l, mol.oop_angles(i,j,k,l));
			}
		}
	}

	printf("torsion!: \n");
	for(int i=0; i < mol.natom; i++)
	{
		for(int j=0; j < i; j++)
		{
			for(int k=0; k < j; k++)
			{
				for(int l=0; l < k; l++)
					if(mol.bond_lengths(i,k) < 4.0 && mol.bond_lengths(k,j) < 4.0 && mol.bond_lengths(k,l) < 4.0)
					printf("%2d-%2d-%2d-%2d %10.6f\n",i,j,k,l, mol.torsion(i,j,k,l));
			}
		}
	}

	double comx = mol.com(0);
	double comy = mol.com(1);
	double comz = mol.com(2);
	printf("Center of mass: \n");
	printf("%8.12f %8.12f %8.12f\n", comx, comy, comz);

	mol.translate(-comx,-comy,-comz);
	printf("Transformed coordinates:\n");
	for(int i = 0; i < mol.natom; i++)
		printf("%2d %10.6f %10.6f %10.6f\n", mol.zvals[i], mol.geom[i][0], mol.geom[i][1], mol.geom[i][2]);

	mol.moi();
	Matrix I(3,3);
	for(int i=0; i < 3; i++)
	{
		for(int j=0; j < 3; j++)
		{
			I(i,j) = mol.moit[i][j];
		}
	}
	std::cout.precision(12);
	std::cout.setf(std::ios::fixed);
	printf("Moment of inertia tensor: \n");
	cout << I << endl;

	//Diagonlize MOI tensor (find its eigenvalues)
	Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
	Matrix evecs = solver.eigenvectors();
	Matrix evals = solver.eigenvalues();
	
	cout << "Principle moments of inertia: " << endl;
	cout << "amu bohr^2: " << endl;
	cout << evals << endl;
	cout << "amu A^2: " << endl;
	cout << evals/(1.8897259886*1.8897259886) << endl;

	//Molecule classification
	if(mol.natom ==2) cout << "\nMolecule is diatomic.\n";
	else if(evals(0) < 1e-4) cout << "\nMolecule is linear.\n";
	else if((fabs(evals(0) - evals(1)) < 1e-4) && (fabs(evals(1) - evals(2)) < 1e-4))
    cout << "\nMolecule is a spherical top.\n";
    else if((fabs(evals(0) - evals(1)) < 1e-4) && (fabs(evals(1) - evals(2)) > 1e-4))
    cout << "\nMolecule is an oblate symmetric top.\n";
    else if((fabs(evals(0) - evals(1)) > 1e-4) && (fabs(evals(1) - evals(2)) < 1e-4))
    cout << "\nMolecule is a prolate symmetric top.\n";
    else cout << "\nMolecule is an asymmetric top.\n";

	
	//Rotational constants
	cout << "Rotational constants: " << endl;
	cout << "MHz: " << endl;
	cout <<	"A = " << mol.rot_const(evals(0))*1e-6 << "  B = " << mol.rot_const(evals(1))*1e-6 << "  C = " << mol.rot_const(evals(2))*1e-6 << endl;
	cout << "cm-1: " << endl;
	cout <<	"A = " << mol.rot_const(evals(0))/2.99792458E10 << "  B = " << mol.rot_const(evals(1))/2.99792458E10 << "  C = " << mol.rot_const(evals(2))/2.99792458E10 << endl;

	hello();
	return 0;
}

