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
	//cout << setprecision(6) << "\n \n" << "mo fock" << "\n";
	//cout << mo_F << "\n \n";

	//Dipole moment
//	cout << "Dipole moments: \n";
//	cout << "x: " << dipole(0) << endl;
//	cout << "y: " << dipole(1) << endl;
//	cout << "z: " << dipole(2) << endl;

}


void hf::mp2_proc()
{
	cout << setprecision(12);
	cout << "~~~MP2~~~" << endl;

//	noddy();
	smarty();
	mp2_energy();

	cout << "Escf = " << e_total << endl;
	cout << "Emp2 = " << Emp2 << endl;
	cout << "Etot = " << e_total + Emp2 << endl;
}

void hf::cc_proc(double ccEconv, double ccRMSconv)
{
	cout << "~~~CCSD~~~" << endl;
	spin_orbit();
	spin_mo_fock();
	init_amp();

	bool cc_rms_converged = 0;
	bool cc_e_converged = 0;
	cc_iterate();
	cout << "iter = " << setw(3) << iter << "  Ecc =    " << cc_e << endl;
	do
	{
		cc_iterate();
		cout << "iter = " << setw(3) << iter << "  Ecc =    " << cc_e << endl;
		if(abs(cc_e - cc_e_old) < ccEconv)
			cc_e_converged = 1;
		cc_rms_converged = cc_conv_check(ccRMSconv);
	}
	while(!cc_rms_converged || !cc_e_converged);

	cout << "CC iterations converged" << endl;
	cout << "Ecc  = "  << cc_e << endl;
	cout << "Etot = "  << e_total + cc_e << endl;

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

bool hf::cc_conv_check(double cc_rms_conv)
{
	double cc_rms = 0;
	for(int p=0; p < nmo; p++)
		for(int q=0; q < nmo; q++)
		{
			cc_rms += (T1[p][q] - T1_old[p][q])*(T1[p][q] - T1_old[p][q]);
			for(int r=0; r < nmo; r++)
				for(int s=0; s < nmo; s++)
				{
					cc_rms += (T2[p][q][r][s] - T2_old[p][q][r][s]) * (T2[p][q][r][s] - T2_old[p][q][r][s]);
				}
		}

	cc_rms = sqrt(cc_rms);
	if(cc_rms < cc_rms_conv)
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

void hf::noddy()
{
	int i, j, k, l, ijkl;
	int p, q, r, s, pq, rs, pqrs;

	for(i=0,ijkl=0; i < orbitals; i++) {
		for(j=0; j <= i; j++) {
			for(k=0; k <= i; k++) {
				for(l=0; l <= (i==k ? j : k); l++,ijkl++) {

					for(p=0; p < orbitals; p++) {
						for(q=0; q < orbitals; q++) {
							pq = INDEX(p,q);
							for(r=0; r < orbitals; r++) {
								for(s=0; s < orbitals; s++) {
									rs = INDEX(r,s);
									pqrs = INDEX (pq,rs);

									TEI_MO[ijkl] += C(p,i) * C(q,j) * C(r,k) * C(s,l) * TEI[pqrs];
								}
							}
						}
					}
				}
			}
		}
	}
}

void hf::smarty()
{
	Matrix X(orbitals,orbitals);
	Matrix Y(orbitals,orbitals);
	Matrix TMP(orbitals*(orbitals+1)/2,orbitals*(orbitals+1)/2);

	//Loop over every k and l pair
	for(int k=0; k < orbitals; k++)
	{
		for(int l=0; l <= k; l++)
		{
			int kl = INDEX(k,l);

			for(int m=0; m < orbitals; m++)
			{
				//Every j and m pair
				for(int j=0; j < orbitals; j++)
				{
					X(m,j) = 0;

					//Transform i to m
					for(int i=0; i < orbitals; i++)
					{
						int ij = INDEX(i,j);
						int ijkl = INDEX(ij,kl);	
						X(m,j) += C(i,m)*TEI[ijkl];
					}	
				}
				
				//Every m and n pair (w/ symmetry)
				for(int n=0; n <= m; n++)
				{
					Y(m,n) = 0;

					//j to n
					for(int j=0; j < orbitals; j++)
					{
						int jm = INDEX(j,m);
						int jmkl = INDEX(jm,kl);
						Y(m,n) += C(j,n)*X(m,j);
					}
					Y(n,m) = Y(m,n);

					//Store half-transformed ints
					int mn = INDEX(m,n);
					TMP(mn,kl) = Y(m,n);
				}
			}
		}
	}

	
	//Loop over every m and n pair
	for(int m=0; m < orbitals; m++)
	{
		for(int n=0; n <= m; n++)
		{
			int mn = INDEX(m,n);


			for(int p=0; p < orbitals; p++)
			{
				//Loop over p and l
				for(int l=0; l < orbitals; l++)
				{
					X(p,l) = 0;
					
					//k to p
					for(int k=0; k < orbitals; k++)
					{
						int kl = INDEX(k,l);
						X(p,l) += C(k,p)*TMP(mn,kl);
					}
				}


				//p and q pairs(w/ sym)
				for(int q=0; q <= p; q++)
				{
					//prep array storage
					int pq = INDEX(p,q);
					int mnpq = INDEX(mn,pq);
					TEI_MO[mnpq] = 0;

					//l to q and store!
					for(int l=0; l < orbitals; l++)
					{
						TEI_MO[mnpq] += C(l,q) * X(p,l);
					}
				}
			}
		}
	}
	

}

void hf::mp2_energy()
{
	Emp2 = 0.0;
	int i, j, a, b, ja, jb, ia, ib, iajb, ibja;
	int ndocc = nelec/2;	

	for(i=0; i < ndocc; i++)
	{
		for(a=ndocc; a < orbitals; a++)
		{
			ia = INDEX(i,a);
			for(j=0; j < ndocc; j++)
			{
				ja = INDEX(j,a);
				for(b=ndocc; b < orbitals; b++)
				{
					jb = INDEX(j,b);
					ib = INDEX(i,b);
					iajb = INDEX(ia,jb);
					ibja = INDEX(ib,ja);
					Emp2 += TEI_MO[iajb] * (2 * TEI_MO[iajb] - TEI_MO[ibja])/(mo_F(i,i) + mo_F(j,j) - mo_F(a,a) - mo_F(b,b));
				}
			}
		}
	}


}

void hf::spin_orbit()
{
	int p, q, r, s, pr, qs, prqs, ps, qr, psqr;
	double value1, value2;
	for(int p=0; p < nmo; p++)
		for(int q=0; q < nmo; q++)
			for(int r=0; r < nmo; r++)
				for(int s=0; s < nmo; s++)
				{
					pr = INDEX(p/2,r/2);
					qs = INDEX(q/2,s/2);
					prqs = INDEX(pr,qs);
					value1 = TEI_MO[prqs] * (p%2 == r%2) * (q%2 == s%2);
					ps = INDEX(p/2,s/2);
					qr = INDEX(q/2,r/2);
					psqr = INDEX(ps,qr);
					value2 = TEI_MO[psqr] * (p%2 == s%2) * (q%2 == r%2);
					ints[p][q][r][s] = value1 - value2;
				}
}

void hf::spin_orbit_CI()
{
	int p, q, r, s, pr, qs, prqs, ps, qr, psqr;
	double value1, value2;
	for(int p=0; p < nmo; p++)
		for(int q=0; q < nmo; q++)
			for(int r=0; r < nmo; r++)
				for(int s=0; s < nmo; s++)
				{
					pr = INDEX(p/2,r/2);
					qs = INDEX(q/2,s/2);
					prqs = INDEX(pr,qs);
					value1 = TEI_MO[prqs] * (p%2 == r%2) * (q%2 == s%2);
					ints[p][q][r][s] = value1;
				}


	//spin_MO fock part
	Matrix mo_H(orbitals,orbitals);
	for(int i=0; i < orbitals; i++)
	{
		for(int j=0; j < orbitals; j++)
		{
			mo_H(i,j) = 0;
			for(int k=0;  k < orbitals; k++)
			{
				for(int l=0; l < orbitals; l++)
				{
					mo_H(i,j) += C(k,j) * C(l,i) * Hcore(k,l);
				}
			}
		}
	}

	for(int p=0; p < nmo; p++)
	{
		for(int q=0; q < nmo; q++)
		{
			smo_F[p][q] = mo_H(p/2,q/2)*(p%2 == q%2);
			for(int m=0; m < nelec; m++)
			{
				smo_F[p][q] += ints[p][m][q][m] - ints[p][m][m][q];
			}
		}
	}
}


void hf::spin_mo_fock()
{
	Matrix mo_H(orbitals,orbitals);
	for(int i=0; i < orbitals; i++)
	{
		for(int j=0; j < orbitals; j++)
		{
			mo_H(i,j) = 0;
			for(int k=0;  k < orbitals; k++)
			{
				for(int l=0; l < orbitals; l++)
				{
					mo_H(i,j) += C(k,j) * C(l,i) * Hcore(k,l);
				}
			}
		}
	}

	for(int p=0; p < nmo; p++)
	{
		for(int q=0; q < nmo; q++)
		{
			smo_F[p][q] = mo_H(p/2,q/2)*(p%2 == q%2);
			for(int m=0; m < nelec; m++)
			{
				smo_F[p][q] += ints[p][m][q][m];
			}
		}
	}

}

void hf::test_mo_fock()
{
	Matrix testfock(nmo,nmo);
	for(int p=0; p < nmo; p++)
	{
		for(int q=0; q < nmo; q++)
		{

			if(p == q)
			{
				smo_F[p][q] = mo_F(p/2,q/2);
			}
			else
			{
				smo_F[p][q] = 0.0;

			}
			
//			smo_F[p][q] = mo_F(p/2,q/2);
			testfock(p,q) = smo_F[p][q];
		}
	}
	cout << testfock << endl;

}


void hf::init_amp()
{
	double emp2_ccsd = 0;
	for(int i=0; i < nelec; i++)
	{
		for(int a=nelec; a < nmo; a++)
		{
			T1[i][a] = 0;
			for(int j=0; j < nelec; j++)
			{
				for(int b=nelec; b < nmo; b++)
				{
					double orb_e = smo_F[i][i] + smo_F[j][j] - smo_F[a][a] - smo_F[b][b];
					T2[i][j][a][b] = ints[i][j][a][b]/orb_e;
					emp2_ccsd += ints[i][j][a][b] * T2[i][j][a][b];
				}
			}
		}
	}

	emp2_ccsd = emp2_ccsd/4.0;

	cout << "Emp2ccsd = " << emp2_ccsd << endl;

	for(int i=0; i < nelec; i++)
	{
		for(int a=nelec; a < nmo; a++)
		{
			Dia[i][a] = smo_F[i][i] - smo_F[a][a];
			for(int j=0; j < nelec; j++)
			{
				for(int b=nelec; b < nmo; b++)
				{
					Dijab[i][j][a][b] = smo_F[i][i] + smo_F[j][j] - smo_F[a][a] - smo_F[b][b];
				}
			}
		}
	}

}

void hf::cc_inter()
{
	//Tau
	for(int i=0; i < nelec; i++)
	{
		for(int a=nelec; a < nmo; a++)
		{
			for(int j=0; j < nelec; j++)
			{
				for(int b=nelec; b < nmo; b++)
				{
					tilde_tau[i][j][a][b] = T2[i][j][a][b] + (T1[i][a]*T1[j][b] - T1[i][b]*T1[j][a])/2.0;
					tau[i][j][a][b] = T2[i][j][a][b] + T1[i][a]*T1[j][b] - T1[i][b]*T1[j][a];

				}
			}
		}
	}

	//Fae
	for(int a=nelec; a < nmo; a++)
	{
		for(int e=nelec; e < nmo; e++)
		{
			Fae[a][e] = (1 - (a == e))*smo_F[a][e];
			for(int m=0; m < nelec; m++)
			{
				Fae[a][e] -= 0.5*smo_F[m][e]*T1[m][a];
				for(int f=nelec; f < nmo; f++)
				{
					Fae[a][e] += T1[m][f] * ints[m][a][f][e];
					for(int n=0; n < nelec; n++)
					{
						Fae[a][e] -= 0.5*tilde_tau[m][n][a][f] * ints[m][n][e][f];
					}
				}
			}
		}
	}

	//Fmi
	for(int m=0; m < nelec; m++)
	{
		for(int i=0; i < nelec; i++)
		{
			Fmi[m][i] = (1 - (m==i))*smo_F[m][i];
			for(int e=nelec; e < nmo; e++)
			{
				Fmi[m][i] += 0.5*T1[i][e]*smo_F[m][e];
				for(int n=0; n < nelec; n++)
				{
					Fmi[m][i] += T1[n][e]*ints[m][n][i][e];
					for(int f=nelec; f < nmo; f++)
					{
						Fmi[m][i] += 0.5*tilde_tau[i][n][e][f]*ints[m][n][e][f];
					}
				}
			}
		}
	}	

	//Fme
	for(int m=0; m < nelec; m++)
	{
		for(int e=nelec; e < nmo; e++)
		{
			Fme[m][e] = smo_F[m][e];
			for(int n=0; n < nelec; n++)
			{
				for(int f=nelec; f < nmo; f++)
				{
					Fme[m][e] += T1[n][f] * ints[m][n][e][f];
				}
			}
		}
	}
	
	//Wmnij
	for(int m=0; m < nelec; m++)
	{
		for(int n=0; n < nelec; n++)
		{
			for(int i=0; i < nelec; i++)
			{
				for(int j=0; j < nelec; j++)
				{
					Wmnij[m][n][i][j] = ints[m][n][i][j];
					for(int e=nelec; e < nmo; e++)
					{
						Wmnij[m][n][i][j] += T1[j][e] * ints[m][n][i][e]; 
						Wmnij[m][n][i][j] -= T1[i][e] * ints[m][n][j][e];
						for(int f=nelec; f < nmo; f++)
						{
							Wmnij[m][n][i][j] += 0.25*tau[i][j][e][f] * ints[m][n][e][f];
						}
					}
				}
			}
		}
	}

	//Wabef
	for(int a=nelec; a < nmo; a++)
	{
		for(int b=nelec; b < nmo; b++)
		{
			for(int e=nelec; e < nmo; e++)
			{
				for(int f=nelec; f < nmo; f++)
				{
					Wabef[a][b][e][f] = ints[a][b][e][f];
					for(int m=0; m < nelec; m++)
					{
						Wabef[a][b][e][f] -= T1[m][b] * ints[a][m][e][f];
						Wabef[a][b][e][f] += T1[m][a] * ints[b][m][e][f];
						for(int n=0; n < nelec; n++)
						{
							Wabef[a][b][e][f] += (tau[m][n][a][b] * ints[m][n][e][f]) / 4.0;
						}
					}
				}	
			}
		}
	}
	//Wmbej
	
	for(int m=0; m < nelec; m++)
	{
		for(int b=nelec; b < nmo; b++)
		{
			for(int e=nelec; e < nmo; e++)
			{
				for(int j=0; j < nelec; j++)
				{
					Wmbej[m][b][e][j] = ints[m][b][e][j];
					for(int f=nelec; f < nmo; f++)
					{
						Wmbej[m][b][e][j] += T1[j][f] * ints[m][b][e][f];
					}
					for(int n=0; n < nelec; n++)
					{
						Wmbej[m][b][e][j] -= T1[n][b] * ints[m][n][e][j];
						for(int f=nelec; f < nmo; f++)
						{
							Wmbej[m][b][e][j] -= (0.5*T2[j][n][f][b] + T1[j][f]*T1[n][b]) * ints[m][n][e][f];
						}
					}
				}
			}
		}
	}

}

void hf::amp_up()
{
	//T1 and T2 old
	for(int p=0; p < nmo; p++)
	{
		for(int q=0; q < nmo; q++)
		{
			T1_old[p][q] = T1[p][q];
			for(int r=0; r < nmo; r++)
			{
				for(int s=0; s < nmo; s++)
				{
					T2_old[p][q][r][s] = T2[p][q][r][s];
				}
			}
		}
	}
	
	//T1
	for(int i=0; i < nelec; i++)
	{
		for(int a=nelec; a < nmo; a++)
		{
			T1[i][a] = smo_F[i][a];
			for(int e=nelec; e < nmo; e++)
			{
				T1[i][a] += T1_old[i][e] * Fae[a][e];
			}
			for(int m=0; m < nelec; m++)
			{
				T1[i][a] -= T1_old[m][a] * Fmi[m][i];
				for(int e=nelec; e < nmo; e++)
				{
					T1[i][a] += T2_old[i][m][a][e] * Fme[m][e];
				}
			}
			for(int n=0; n < nelec; n++)
			{
				for(int f=nelec; f < nmo; f++)
				{
					T1[i][a] -= T1_old[n][f] * ints[n][a][i][f];
				}
			}
			for(int m=0; m < nelec; m++)
			{
				for(int e=nelec; e < nmo; e++)
				{
					for(int f=nelec; f < nmo; f++)
					{
						T1[i][a] -= 0.5*T2_old[i][m][e][f] * ints[m][a][e][f];
					}
					for(int n=0; n < nelec; n++)
					{
						T1[i][a] -= 0.5*T2_old[m][n][a][e] * ints[n][m][e][i];
					}
				}
			}
			T1[i][a] /= Dia[i][a];	

			//T2
			for(int j=0; j < nelec; j++)
			{
				for(int b=nelec; b < nmo; b++)
				{
					T2[i][j][a][b] = ints[i][j][a][b];
					for(int e=nelec; e < nmo; e++)
					{
						//hail mary
						T2[i][j][a][b] += T2_old[i][j][a][e]*Fae[b][e];
						T2[i][j][a][b] -= T2_old[i][j][b][e]*Fae[a][e];
						for(int m=0; m < nelec; m++)
						{
							T2[i][j][a][b] -= T2_old[i][j][a][e] * 0.5 * T1_old[m][b] * Fme[m][e];
							T2[i][j][a][b] += T2_old[i][j][b][e] * 0.5 * T1_old[m][a] * Fme[m][e];

						}
					}
					for(int m=0; m < nelec; m++)
					{
						T2[i][j][a][b] -= T2_old[i][m][a][b] * Fmi[m][j];
						T2[i][j][a][b] += T2_old[j][m][a][b] * Fmi[m][i];
						for(int e=nelec; e < nmo; e++)
						{
							//*crosses fingers*
							T2[i][j][a][b] -= 0.5 * T2_old[i][m][a][b] * T1_old[j][e] * Fme[m][e];
							T2[i][j][a][b] += 0.5 * T2_old[j][m][a][b] * T1_old[i][e] * Fme[m][e];
						}

						for(int n=0; n < nelec; n++)
						{
							T2[i][j][a][b] += 0.5 * tau[m][n][a][b] * Wmnij[m][n][i][j];
						}
					}
					
					for(int e=nelec; e < nmo; e++)
					{
						for(int f=nelec; f < nmo; f++)
						{
							T2[i][j][a][b] += 0.5 * tau[i][j][e][f] * Wabef[a][b][e][f];
						}
					}

					//Last two terms you can do it Megan!!!

					for(int m=0; m < nelec; m++)
					{
						for(int e=nelec; e < nmo; e++)
						{
							//let's test it here!
							T2[i][j][a][b] += T2_old[i][m][a][e] * Wmbej[m][b][e][j] - (T1_old[i][e] * T1_old[m][a] * ints[m][b][e][j]);


							T2[i][j][a][b] -= T2_old[i][m][b][e] * Wmbej[m][a][e][j] - (T1_old[i][e] * T1_old[m][b] * ints[m][a][e][j]);


							T2[i][j][a][b] -= T2_old[j][m][a][e] * Wmbej[m][b][e][i] - (T1_old[j][e] * T1_old[m][a] * ints[m][b][e][i]);


							T2[i][j][a][b] += T2_old[j][m][b][e] * Wmbej[m][a][e][i] - (T1_old[j][e] * T1_old[m][b] * ints[m][a][e][i]);

						}
					}

					//Last two for real term!

					for(int e=nelec; e < nmo; e++)
					{
						T2[i][j][a][b] += T1_old[i][e] * ints[a][b][e][j];
						T2[i][j][a][b] -= T1_old[j][e] * ints[a][b][e][i];
					}

					for(int m=0; m < nelec; m++)
					{
						T2[i][j][a][b] -= T1_old[m][a] * ints[m][b][i][j];
						T2[i][j][a][b] += T1_old[m][b] * ints[m][a][i][j];
					}

					T2[i][j][a][b] /= Dijab[i][j][a][b];
				}
			}
		}
	}

}

void hf::cc_energy()
{
	cc_e_old = cc_e;
	cc_e = 0;

	for(int i=0; i < nelec; i++)
	{
		for(int a=nelec; a < nmo; a++)
		{
			cc_e += smo_F[i][a] * T1[i][a];
			for(int j=0; j < nelec; j++)
			{
				for(int b=nelec; b < nmo; b++)
				{
					cc_e += 0.25*ints[i][j][a][b] * T2[i][j][a][b];
					cc_e += 0.5*ints[i][j][a][b] * T1[i][a] * T1[j][b];
				}
			}
		}
	}



}

void hf::cc_iterate()
{
	cc_inter();
	amp_up();
	cc_energy();
	iter += 1;
}


void hf::full_pt()
{
	double Dijkabc;
	double T3D;
	double T3C;
	//Dijkabc
	e_pt = 0;
	for(int i=0; i < nelec; i++)
		for(int j=0; j < nelec; j++)
			for(int k=0; k < nelec; k++)
				for(int a=nelec; a < nmo; a++)
					for(int b=nelec; b < nmo; b++)
						for(int c=nelec; c < nmo; c++)
						{
							//Dijkabc
							Dijkabc = smo_F[i][i] + smo_F[j][j] + smo_F[k][k] - smo_F[a][a] - smo_F[b][b] - smo_F[c][c];

							//Disconnected triples
							//ijk permutations first
							T3D = T1[i][a] * ints[j][k][b][c];
							T3D -= T1[i][b] * ints[j][k][a][c];
							T3D -= T1[i][c] * ints[j][k][b][a];

							//jik permutations
							T3D -= T1[j][a] * ints[i][k][b][c];
							T3D += T1[j][b] * ints[i][k][a][c];
							T3D += T1[j][c] * ints[i][k][b][a];

							//kji permutations
							T3D -= T1[k][a] * ints[j][i][b][c];
							T3D += T1[k][b] * ints[j][i][a][c];
							T3D += T1[k][c] * ints[j][i][b][a];

							T3D /= Dijkabc;

							//Connected triples
							T3C = 0;
							for(int e=nelec; e < nmo; e++)
							{
								//ijk
								T3C += T2[j][k][a][e] * ints[e][i][b][c];
								T3C -= T2[j][k][b][e] * ints[e][i][a][c];
								T3C -= T2[j][k][c][e] * ints[e][i][b][a];

								//jik
								T3C -= T2[i][k][a][e] * ints[e][j][b][c];
								T3C += T2[i][k][b][e] * ints[e][j][a][c];
								T3C += T2[i][k][c][e] * ints[e][j][b][a];

								//kji
								T3C -= T2[j][i][a][e] * ints[e][k][b][c];
								T3C += T2[j][i][b][e] * ints[e][k][a][c];
								T3C += T2[j][i][c][e] * ints[e][k][b][a];
							}
							for(int m=0; m < nelec; m++)
							{
								//ijk
								T3C -= T2[i][m][b][c] * ints[m][a][j][k];
								T3C += T2[i][m][a][c] * ints[m][b][j][k];
								T3C += T2[i][m][b][a] * ints[m][c][j][k];

								//jik
								T3C += T2[j][m][b][c] * ints[m][a][i][k];
								T3C -= T2[j][m][a][c] * ints[m][b][i][k];
								T3C -= T2[j][m][b][a] * ints[m][c][i][k];

								//kji
								T3C += T2[k][m][b][c] * ints[m][a][j][i];
								T3C -= T2[k][m][a][c] * ints[m][b][j][i];
								T3C -= T2[k][m][b][a] * ints[m][c][j][i];

							}
							T3C /= Dijkabc;

							e_pt += (T3C*Dijkabc*(T3C + T3D))/36.0;
						}

	cout << "~~~(T)~~~" << endl;
	cout << "E(T) = " << e_pt << endl;
	cout << "ECCSD(T) = " << e_total + cc_e + e_pt << endl;

}

void hf::CI()
{
	//Subtract nelec from a to get relative index of virtual orbital
	int count = 0;
	for(int i=0; i < nelec; i++)
		for(int a=nelec; a < nmo; a++)
		{
//			int ia = INDEX(i,a - nelec);
			int ia = count;
			count++;
			int countb = 0;
			for(int j=0; j < nelec; j++)
				for(int b=nelec; b < nmo; b++)
				{
					int jb = countb;
					HCI(ia,jb) = smo_F[a][b]*(i==j) - smo_F[i][j]*(a==b) + ints[a][j][i][b];
					countb++;
				}
		}
//	cout << HCI << endl;


	//Diagonalize
	Eigen::SelfAdjointEigenSolver<Matrix> solver(HCI);
	Matrix CIS_XS = solver.eigenvectors();
	Matrix CIS_E = solver.eigenvalues();

	//Pretty print
	cout << "~~CIS~~" << endl;
	cout << "Excitation energies:" << endl;
	cout << CIS_E << endl;
}


void hf::CI_spin()
{

	int count = 0;
	for(int i=0; i < nelec/2; i++)
		for(int a=nelec/2; a < orbitals; a++)
		{
			int ia = count;	
			count++;
			int countb = 0;
			for(int j=0; j < nelec/2; j++)
				for(int b=nelec/2; b < orbitals; b++)
				{
					int jb = countb;
	//				cout << ia << " " << jb << endl;
					HCI_spin(ia,jb) = mo_F(a,b)*(i==j) - mo_F(i,j)*(a==b) + 2*ints[a*2][j*2][i*2][b*2] - ints[a*2][j*2][b*2][i*2]; 
					countb++;
				}
		}


//	cout << HCI_spin << endl;
	//diagonlize??
	Eigen::SelfAdjointEigenSolver<Matrix> solver(HCI_spin);
	Matrix CIS_XS = solver.eigenvectors();
	Matrix CIS_E = solver.eigenvalues();


	cout << "\n~~Spin-Adapted CIS~~"<< endl;
	cout << "Excitation energies (Hartree):" << endl;
	cout << "Singlets:" << endl;
	cout << CIS_E << endl;

	count = 0;
	for(int i=0; i < nelec/2; i++)
		for(int a=nelec/2; a < orbitals; a++)
		{
			int ia = count;	
			count++;
			int countb = 0;
			for(int j=0; j < nelec/2; j++)
				for(int b=nelec/2; b < orbitals; b++)
				{
					int jb = countb;
	//				cout << ia << " " << jb << endl;
					HCI_spin(ia,jb) = mo_F(a,b)*(i==j) - mo_F(i,j)*(a==b) - ints[a*2][j*2][b*2][i*2]; 
					countb++;
				}
		}
	Eigen::SelfAdjointEigenSolver<Matrix> solver2(HCI_spin);
	 CIS_XS = solver2.eigenvectors();
	 CIS_E = solver2.eigenvalues();
	cout << "Triplets:" << endl;
	cout << CIS_E << endl;

}


void hf::RPA()
{

	//Calculate how much to adjust position for B and A-
	int vert = nmo - nelec;
	int vn = vert*nelec;

	int	count = 0;
	for(int i=0; i < nelec; i++)
		for(int a=nelec; a < nmo; a++)
		{
			int ia = count;
			count++;
			int countb = 0;
			for(int j=0; j < nelec; j++)
				for(int b=nelec; b < nmo; b++)
				{
					int jb = countb;
					//A matrix
					//So <aj||ib> = <aj|ib> - <ab|ji>
					HRPA(ia,jb) = smo_F[a][b]*(i==j) - smo_F[i][j]*(a==b) + ints[a][j][i][b] - ints[a][j][b][i];
					//-A
					HRPA(ia+vn,jb+vn) = -HRPA(ia,jb);
					//B
					HRPA(ia,jb+vn) = -ints[a][b][i][j] + ints[a][b][j][i];
					//-B
					HRPA(ia+vn,jb) = -HRPA(ia,jb+vn);
					countb++;
				}
		}

	//Figure out printing later
	cout << "\n~~TDHF/RPA~~" << endl;
	cout << HRPA.eigenvalues().real() << endl;

}

void hf::RPA_reduced()
{
	int vert = nmo - nelec;
	int vn = vert*nelec;
	Matrix A(vn,vn);
	Matrix B(vn,vn);
	int	count = 0;
	for(int i=0; i < nelec; i++)
		for(int a=nelec; a < nmo; a++)
		{
			int ia = count;
			count++;
			int countb = 0;
			for(int j=0; j < nelec; j++)
				for(int b=nelec; b < nmo; b++)
				{
					int jb = countb;
					A(ia,jb) = smo_F[a][b]*(i==j) - smo_F[i][j]*(a==b) + ints[a][j][i][b] - ints[a][j][b][i];
					B(ia,jb) = ints[a][b][i][j] - ints[a][b][j][i];
					countb++;
				}
		}
	HRPA_reduced = (A + B)*(A-B);
	cout << "\n~~TDHF/RPA~~" << endl;
	Matrix rip = HRPA_reduced.eigenvalues().real().array().sqrt();
	vector<double> v(vn);
//	Matrix rip = HRPA_reduced.eigenvalues().real().array().sqrt();
	for(int i=0; i < vn; i++)
		v[i] = rip(i);
	sort(v.begin(), v.end());
	cout << "Excitation energies (Hartree)" << endl;
	for(int i=0; i < vn; i++)
		cout << v[i] << endl;
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
	TEI_MO = new double[ijkl_max+1];

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


	//Allocate memory for CCSD stuff just doing it here because why wouldn't I
	nmo = orbitals*2;
	ints = new double*** [nmo];
	smo_F = new double* [nmo];
	T2 = new double*** [nmo];
	T2_old = new double*** [nmo];
	T1 = new double* [nmo];
	T1_old = new double* [nmo];
	Fae = new double* [nmo];
	Fmi = new double* [nmo];
	Fme = new double* [nmo];
	Wmnij = new double*** [nmo];
	Wabef = new double*** [nmo];
	Wmbej = new double*** [nmo];
	tilde_tau = new double*** [nmo];
	tau = new double*** [nmo];
	Dia = new double* [nmo];
	Dijab = new double*** [nmo];

	//(T) stuff
//	Dijkabc = new double***** [nmo];
//	T3D = new double***** [nmo];
//	T3C = new double***** [nmo];


	for(int i=0; i < nmo; i++)
	{
		ints[i] = new double** [nmo];
		T2[i] = new double** [nmo];
		T2_old[i] = new double** [nmo];
		smo_F[i] = new double[nmo];
		T1[i] = new double[nmo];
		T1_old[i] = new double[nmo];
 		Fae[i] = new double[nmo];
 		Fmi[i] = new double[nmo];
 		Fme[i] = new double[nmo];
 		Wmnij[i] = new double** [nmo];
 		Wabef[i] = new double** [nmo];
 		Wmbej[i] = new double** [nmo];
		tilde_tau[i] = new double** [nmo];
		tau[i] = new double** [nmo];
		Dia[i] = new double[nmo];
		Dijab[i] = new double** [nmo];
//		Dijkabc[i] = new double**** [nmo];
//		T3D[i] = new double**** [nmo];
//		T3C[i] = new double**** [nmo];
		for(int j=0; j < nmo; j++)
		{
			ints[i][j] = new double* [nmo];	
			T2[i][j] = new double* [nmo];	
			T2_old[i][j] = new double* [nmo];	
 			Wmnij[i][j] = new double* [nmo];
 			Wabef[i][j] = new double* [nmo];
 			Wmbej[i][j] = new double* [nmo];
			tilde_tau[i][j] = new double* [nmo];
			tau[i][j] = new double* [nmo];
			Dijab[i][j] = new double* [nmo];
//			Dijkabc[i][j]  = new double*** [nmo];
//			T3D[i][j]  = new double*** [nmo];
//			T3C[i][j]  = new double*** [nmo];
			for(int k=0; k < nmo; k++)
			{
				ints[i][j][k] = new double[nmo];
				T2[i][j][k] = new double[nmo];
				T2_old[i][j][k] = new double[nmo];
 				Wmnij[i][j][k] = new double[nmo];
 				Wabef[i][j][k] = new double[nmo];
 				Wmbej[i][j][k] = new double[nmo];
				tilde_tau[i][j][k] = new double[nmo];
				tau[i][j][k] = new double[nmo];
				Dijab[i][j][k] = new double[nmo];
//				Dijkabc[i][j][k] = new double** [nmo];
//				T3D[i][j][k] = new double** [nmo];
//				T3C[i][j][k] = new double** [nmo];
				/*
				for(int l=0; l < nmo; l++)
				{
					Dijkabc[i][j][k][l] = new double* [nmo];
					T3D[i][j][k][l]  = new double* [nmo];
					T3C[i][j][k][l]  = new double* [nmo];
					for(int m=0; m < nmo; m++)
					{
						Dijkabc[i][j][k][l][m]  = new double[nmo];
						T3D[i][j][k][l][m] = new double[nmo];
						T3C[i][j][k][l][m] = new double[nmo];
					}
				}
				*/
			}
		}
	}

	cc_e = 0;
	iter = 0;
	e_pt = 0;

	//Initialize HCI
	int vert = nmo - nelec;
	int vn = vert*nelec;
	int vert_spat = orbitals - nelec/2; 
	int vn_spat = vert_spat * (nelec/2);	
	HCI.resize(vn,vn);
	HCI_spin.resize(vn_spat,vn_spat);

	//RPA
	HRPA.resize(vn*2,vn*2);
	HRPA_reduced.resize(vn,vn);
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

	for(int i=0; i < nmo; i++)
	{
		for(int j=0; j < nmo; j++)
		{
			for(int k=0; k < nmo; k++)
			{
				delete[] ints[i][j][k];
				delete[] T2[i][j][k];
				delete[] T2_old[i][j][k];
				delete[] Wmnij[i][j][k];
				delete[] Wabef[i][j][k];
				delete[] Wmbej[i][j][k];
				delete[] tilde_tau[i][j][k];
				delete[] tau[i][j][k];
				delete[] Dijab[i][j][k];
				for(int l=0; l < nmo; l++)
				{
					for(int m=0; m < nmo; m++)
					{
//						delete[] Dijkabc[i][j][k][l][m];
//						delete[] T3D[i][j][k][l][m];
//						delete[] T3C[i][j][k][l][m];
					}
//					delete[] Dijkabc[i][j][k][l];
//					delete[] T3D[i][j][k][l];
//					delete[] T3C[i][j][k][l];
				}
//				delete[] Dijkabc[i][j][k];
//				delete[] T3D[i][j][k];
//				delete[] T3C[i][j][k];
			}
			delete[] ints[i][j];
			delete[] T2[i][j];
			delete[] T2_old[i][j];
			delete[] Wmnij[i][j];
			delete[] Wabef[i][j];
			delete[] Wmbej[i][j];
			delete[] tilde_tau[i][j];
			delete[] tau[i][j];
			delete[] Dijab[i][j];
//			delete[] Dijkabc[i][j];
//			delete[] T3D[i][j];
//			delete[] T3C[i][j];
		}
		delete[] ints[i];
		delete[] smo_F[i];
		delete[] T1[i];
		delete[] T1_old[i];
		delete[] T2[i];
		delete[] T2_old[i];
		delete[] Fae[i];
		delete[] Fmi[i];
		delete[] Fme[i];
		delete[] Wmnij[i];
		delete[] Wabef[i];
		delete[] Wmbej[i];
		delete[] tilde_tau[i];
		delete[] tau[i];
		delete[] Dia[i];
		delete[] Dijab[i];
//		delete[] Dijkabc[i];
//		delete[] T3D[i];
//		delete[] T3C[i];
	}
	delete[] ints;
	delete[] smo_F;
	delete[] T1;
	delete[] T1_old;
	delete[] T2;
	delete[] T2_old;
	delete[] Fae;
	delete[] Fmi;
	delete[] Fme;
	delete[] Wmnij;
	delete[] Wabef;
	delete[] Wmbej;
	delete[] tilde_tau;
	delete[] tau;
	delete[] Dia;
	delete[] Dijab;
//	delete[] Dijkabc;
//	delete[] T3D;
//	delete[] T3C;
}
