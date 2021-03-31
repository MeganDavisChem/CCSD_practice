#include <fstream>
#include <cmath>
#include <cassert>
#include "hf.h"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>

using namespace std;

void hf::fock_procedure()
{
	//Print enuc and Hcore in a nice table
	cout << fixed;

	cout << "enuc" << "\n";
	cout << enuc << "\n \n";
	cout << "Hcore" << "\n";
	cout << Hcore << "\n \n";

	//Form orthogonal transformation (?) matrix and print
	ortho();
	cout << "S-1/2" << "\n";
	cout << S_12 << "\n \n";
	
	//Convergence loop structure
	//Header
	cout << setprecision(12);
	cout << setw(5) << left << "Iter" << '\t' << setw(18) << "E(elec)" << '\t' << setw(18) << "E(tot)" << '\t' << setw(18) << "Delta(E)" << '\t' <<  setw(18) << "RMS(D)" << '\t' << endl;

	int count = 0;
	bool rms_converged = false;
	bool e_converged = false;

	//run loop one time, check convergence after setting e prev  and e elec, use boolean values so it always runs in the right part of the loop
	do
	{
		e_elec_prev = e_elec;
		e_total_prev = e_total;

		//The actual main HF stuff
		guess_density();
		scf_energy();
		new_fock();

		//Output and conv check
		cout << setw(5) << left << count << setw(18) << e_elec << "\t" << setw(18) << e_total << "\t";
		if(count > 0)
		{
			e_converged = e_conv_check(e_conv_crit);
			rms_converged = rms_conv_check(rms_conv_crit);
			cout << setw(18) << deltaE << "\t" << setw(18) << rms << "\t";
		}
		cout << endl;
		count++;
	}
	while(!rms_converged || !e_converged);
	cout << "Converged!" << endl;

	//Now do the Fock in new basis or whatever thing
	mo_fock();
	cout << setprecision(6) << "\n \n" << "mo fock" << "\n";
	cout << mo_F << "\n \n";

	//Dipole moment
//	cout << "Dipole moments: \n";
//	cout << "x: " << dipole(0) << endl;
//	cout << "y: " << dipole(1) << endl;
//	cout << "z: " << dipole(2) << endl;

}

void hf::ortho()
{
	//Convert S to Matrix type
	Matrix Sm(orbitals,orbitals);
	for(int i=0; i < orbitals; i++)
	{
		for(int j=0; j < orbitals; j++)
			Sm(i,j) = S[i][j];
	}

	//Get eigenvalues and vectors
	Eigen::SelfAdjointEigenSolver<Matrix> solver(Sm);
	Matrix Ls = solver.eigenvectors();
	Matrix lambda = solver.eigenvalues();

	//Form pieces for orthoganlization matrix
	Matrix lambda_12 = lambda.array().inverse().sqrt();
	Matrix Ls_trans = Ls.transpose();

	//Actually do the multiplication
	Matrix lambda_12_diag(orbitals,orbitals);
	for(int i=0; i < orbitals; i++)
	{
		for(int j=0; j < orbitals; j++)
		{
			if(i==j)
				lambda_12_diag(i,j) = lambda_12(i);
			else	
				lambda_12_diag(i,j) = 0;
		}
	}
	S_12 = Ls * lambda_12_diag * Ls_trans;
	
	cout << "Ls\n" << Ls << endl;
	cout << "lambda_12\n" << lambda_12 << endl;
    cout << "lamda_12_diag\n" << lambda_12_diag << endl;
    cout << "Ls_trans\n" << Ls_trans << endl;	
	F = Hcore;
}

void hf::guess_density()
{
	D_old = Density;
	F_prime = S_12.transpose() * F * S_12;

	//Get eigenvectors
	Eigen::SelfAdjointEigenSolver<Matrix> solver(F_prime);
	Matrix C_prime = solver.eigenvectors();

	//Transform into original basis
	C = S_12 * C_prime;

	//Density matrix part
	Density.resize(orbitals,orbitals);
	for(int i=0; i < orbitals; i++)
	{
		for(int j=0; j < orbitals; j++)
		{
			Density(i,j) = 0;
			for(int k=0; k < nelec/2; k++)
			{
				Density(i,j) += C(i,k) * C(j,k);
			}
		}
	}
	cout << "C" << C << endl;
}


void hf::scf_energy()
{
	e_elec = 0;	
	for(int i=0; i < orbitals; i++)
	{
		for(int j=0; j < orbitals; j++)
		{
			double temp = Hcore(i,j) + F(i,j);
			e_elec += Density(i,j) * temp;
		}
	}
	e_total = e_elec + enuc;
}

void hf::new_fock()
{
	for(int i=0; i < orbitals; i++)
	{
		for(int j=0; j < orbitals; j++)
		{
			F(i,j) = Hcore(i,j);
			for(int k=0; k < orbitals; k++)
			{
				for(int l=0; l < orbitals; l++)
				{
					int ij = INDEX(i,j);
					int kl = INDEX(k,l);
					int ijkl = INDEX(ij,kl);
					int ik = INDEX(i,k);
					int jl = INDEX(j,l);
					int ikjl = INDEX(ik,jl);

					F(i,j) += Density(k,l) * (2.0 * TEI[ijkl] - TEI[ikjl]);
				}
			}
		}
	}


}

bool hf::e_conv_check(double e_conv)
{
	deltaE = e_elec - e_elec_prev;
	if(abs(deltaE) < e_conv)
	{
		return true;
	}
	else
		return false;
	
}

bool hf::rms_conv_check(double rms_conv)
{
	double t_rms = 0;
	for(int i=0; i < orbitals; i++)
	{
		for(int j=0; j < orbitals; j++)
		{
			double temp = Density(i,j) - D_old(i,j);
			t_rms += temp*temp;
		}

	}
	rms = sqrt(t_rms);

	if(rms < rms_conv)
	{
		return true;
	}
	else
		return false;
}


void hf::mo_fock()
{
	mo_F.resize(orbitals,orbitals);
	for(int i=0; i < orbitals; i++)
	{
		for(int j=0; j < orbitals; j++)
		{
			mo_F(i,j) = 0;
			for(int k=0;  k < orbitals; k++)
			{
				for(int l=0; l < orbitals; l++)
				{
					mo_F(i,j) += C(k,j) * C(l,i) * F(k,l);
				}
			}
		}
	}

}

double hf::dipole(int cart)
{
	cout << cart << endl;
	double mu = 0;
	if(cart == 0)
	{
		for(int i=0; i < orbitals; i++)
		{
			for(int j=0; j < orbitals; j++)
			{
				mu += 2*Density(i,j)*dip_x[i][j];
			}
		}

	}
	else if(cart == 1)
	{
		for(int i=0; i < orbitals; i++)
		{
			for(int j=0; j < orbitals; j++)
			{
				mu += 2*Density(i,j)*dip_y[i][j];
			}
		}
	}
	else if(cart == 2)
	{
		for(int i=0; i < orbitals; i++)
		{
			for(int j=0; j < orbitals; j++)
			{
				mu += 2*Density(i,j)*dip_z[i][j];
			}
		}
	}
	else
	{
		cout << "Sorry!" << endl;
	}

	return mu;
}

hf::hf(string mol, string basis, int elec, double e_conv_crit_input, double rms_conv_crit_input)
{
	e_conv_crit = e_conv_crit_input;
	rms_conv_crit = rms_conv_crit_input;
	e_elec_prev = 0;
	e_total_prev = 0;
	e_elec = 1000;
	e_total = 0;

	nelec = elec;
	string dir = "input/" + mol + "/" + basis + "/";

	//Read enuc
	ifstream is_enuc(dir + "enuc.dat");
	assert(is_enuc.good());
	is_enuc >> enuc;
	is_enuc.close();

	//Read s
	ifstream is_s(dir + "s.dat");
	assert(is_s.good());
	vector<double> s_vec;
	for(double i; is_s >> i;)
		s_vec.push_back(i);
	is_s.close();

	//Allocate memory for S
	int s_size = s_vec[s_vec.size()-2];
	S = new double* [s_size];
	for(int i=0; i < s_size; i++)
		S[i] = new double[s_size];

	//Store S in appropriately constructed matrix...
	for(int i=0; i < s_vec.size()/3; i++)
	{
		int k = s_vec[i*3] - 1;
		int l = s_vec[i*3+1] - 1;
		double x = s_vec[i*3+2];
		S[k][l] = x;
		S[l][k] = x;
	}

	//Read T
	ifstream is_t(dir + "t.dat");
	assert(is_t.good());
	vector<double> t_vec;
	for(double i; is_t >> i;)
		t_vec.push_back(i);
	is_t.close();
	int t_size = t_vec[t_vec.size()-2];
	T = new double* [t_size];
	for(int i=0; i < t_size; i++)
		T[i] = new double[t_size];
	for(int i=0; i < t_vec.size()/3; i++)
	{
		int k = t_vec[i*3] - 1;
		int l = t_vec[i*3+1] - 1;
		double x = t_vec[i*3+2];
		T[k][l] = x;
		T[l][k] = x;
	}

	//Read V
	ifstream is_v(dir + "v.dat");
	assert(is_v.good());
	vector<double> v_vec;
	for(double i; is_v >> i;)
		v_vec.push_back(i);
	is_v.close();
	int v_size = v_vec[v_vec.size()-2];
	V = new double* [v_size];
	for(int i=0; i < v_size; i++)
		V[i] = new double[v_size];
	for(int i=0; i < v_vec.size()/3; i++)
	{
		int k = v_vec[i*3] - 1;
		int l = v_vec[i*3+1] - 1;
		double x = v_vec[i*3+2];
		V[k][l] = x;
		V[l][k] = x;
	}
	//Form Hcore
/*	
	Hcore = new double* [v_size];
	for(int i=0; i < v_size; i++)
	{
		Hcore[i] = new double[v_size];
		for(int j=0; j < v_size; j++)
			Hcore[i][j] = T[i][j] + V[i][j];
	}
	*/
	Hcore.resize(v_size, v_size);
	for(int i=0; i < v_size; i++)
	{
		for(int j=0; j < v_size; j++)
		{
			Hcore(i,j) = T[i][j] + V[i][j];
		}
	}

	orbitals = v_size;

	//Lookup array
	int ij_max = orbitals*(orbitals+1)/2.0 + orbitals;
	int ijkl_max = ij_max*(ij_max+1)/2.0 + ij_max;	
	ioff = new int[ijkl_max];

	ioff[0] = 0;
	for(int i=1; i < ijkl_max; i++)
		ioff[i] = ioff[i-1] + i;


	//TEI mem allocation
	TEI = new double[ijkl_max+1];	

	//open TEI file
	ifstream is_tei(dir + "eri.dat");
	assert(is_tei.good());
	double val;
	for(int i, j, k, l; is_tei >> i >> j >> k >> l >> val;)
	{
		int ij = INDEX(i-1,j-1);
		int kl = INDEX(k-1,l-1);
		int ijkl = INDEX(ij,kl);
		TEI[ijkl] = val;
	}

	D_old.resize(orbitals,orbitals);

	//Allocate memory for dipoles
	dip_x = new double* [orbitals];
	for(int i=0; i < orbitals; i++)
		dip_x[i] = new double[orbitals];

	dip_y = new double* [orbitals];
	for(int i=0; i < orbitals; i++)
		dip_y[i] = new double[orbitals];

	dip_z = new double* [orbitals];
	for(int i=0; i < orbitals; i++)
		dip_z[i] = new double[orbitals];

	//Read in dipole stuff
	ifstream is_x(dir + "mux.dat");
	assert(is_x.good());
	for(int i, j; is_x >> i >> j >> val;)
	{
		dip_x[i-1][j-1] = val;
		dip_x[j-1][i-1] = val;
	}
	is_x.close();

	ifstream is_y(dir + "muy.dat");
	assert(is_y.good());
	for(int i, j; is_y >> i >> j >> val;)
	{
		dip_y[i-1][j-1] = val;
		dip_y[j-1][i-1] = val;
	}
	is_y.close();

	ifstream is_z(dir + "muz.dat");
	assert(is_z.good());
	for(int i, j; is_z >> i >> j >> val;)
	{
		dip_z[i-1][j-1] = val;
		dip_z[j-1][i-1] = val;
	}
	is_z.close();

}

hf::~hf()
{
	//Delete V and T
	for(int i=0; i < orbitals; i++)
	{
		delete[] V[i];
		delete[] T[i];
		delete[] dip_x[i];
		delete[] dip_y[i];
		delete[] dip_z[i];
	}
	delete[] V;
	delete[] T;
	delete[] TEI;
	delete[] ioff;
	delete[] dip_x;
	delete[] dip_y;
	delete[] dip_z;
}
