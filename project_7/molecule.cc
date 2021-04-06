#include <fstream>
#include <cmath>
#include <cassert>
#include "molecule.h"
#include "masses.h"

using namespace std;

void Molecule::print_geom()
{
	for(int i=0; i < natom; i++)
		printf("%d %8.5f %8.5f %8.5f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
}

double Molecule::bond_lengths(int atom1, int atom2)
{
	double bond = sqrt(pow(geom[atom1][0] - geom[atom2][0],2) + pow(geom[atom1][1] - geom[atom2][1],2) + pow(geom[atom1][2] - geom[atom2][2],2)); 
	return bond;
}

void Molecule::mass_weight()
{
	for(int i=0; i < natom*3; i++)
	{
		for(int j=0; j < natom*3; j++)
		H[i][j] = H[i][j]/sqrt(amasses[zvals[i/natom]]*amasses[zvals[j/natom]]);
	}
}

double Molecule::bond_angles(int atom1, int atom2, int atom3)
{
	double angle = acos(unit_vectors(atom2, atom1, 0)*unit_vectors(atom2, atom3, 0) + unit_vectors(atom2, atom1, 1)*unit_vectors(atom2, atom3, 1) + unit_vectors(atom2, atom1, 2)*unit_vectors(atom2, atom3, 2))*180/acos(-1); 
	return angle;
}

double Molecule::unit_vectors(int atom1, int atom2, int cart)
{
	double uv = -(geom[atom2][cart] - geom[atom1][cart])/bond_lengths(atom2,atom1);
	return uv;
}

double Molecule::oop_angles(int atom1, int atom2, int atom3, int atom4)
{
	double uvkj[3];
	double uvkl[3];
	double uvki[3];
	double cross[3];
	for(int i=0; i < 3; i++)
	{
		uvkj[i] = unit_vectors(atom2, atom3, i);
		uvkl[i] = unit_vectors(atom4, atom3, i);
		uvki[i] = unit_vectors(atom1, atom3, i);
	}
	double angle = bond_angles(atom2, atom3, atom4);
	cross[0] = uvkj[1]*uvkl[2] - uvkj[2]*uvkl[1];
	cross[1] = uvkj[2]*uvkl[0] - uvkj[0]*uvkl[2];
	cross[2] = uvkj[0]*uvkl[1] - uvkj[1]*uvkl[0];
	double hmm = (cross[0]*uvki[0] + cross[1]*uvki[1] + cross[2]*uvki[2])/sin(bond_angles(atom2,atom3,atom4)*acos(-1)/180);
	double oop;
	if(hmm < -1)
		oop = asin(-1)*180/acos(-1);
	else if(hmm > 1)
		oop = asin(1)*180/acos(-1);
	else
		oop = asin((cross[0]*uvki[0] + cross[1]*uvki[1] + cross[2]*uvki[2])/sin(bond_angles(atom2,atom3,atom4)*acos(-1)/180))*180/acos(-1);
	return oop;
}

double Molecule::torsion(int atom1, int atom2, int atom3, int atom4)
{
	double uvij[3];
	double uvjk[3];
	double uvkl[3];
	double cross_ij_jk[3];
	double cross_jk_kl[3];
	double ang_ijk = bond_angles(atom1, atom2, atom3)*acos(-1)/180;
	double ang_jkl = bond_angles(atom2, atom3, atom4)*acos(-1)/180;
	for(int i=0; i < 3; i++)
	{
		uvij[i] = unit_vectors(atom2, atom1, i);
		uvjk[i] = unit_vectors(atom3, atom2, i);
		uvkl[i] = unit_vectors(atom4, atom3, i);
	}
	cross_ij_jk[0] = uvij[1]*uvjk[2] - uvij[2]*uvjk[1];
	cross_ij_jk[1] = uvij[2]*uvjk[0] - uvij[0]*uvjk[2];
	cross_ij_jk[2] = uvij[0]*uvjk[1] - uvij[1]*uvjk[0];

	cross_jk_kl[0] = uvjk[1]*uvkl[2] - uvjk[2]*uvkl[1];
	cross_jk_kl[1] = uvjk[2]*uvkl[0] - uvjk[0]*uvkl[2];
	cross_jk_kl[2] = uvjk[0]*uvkl[1] - uvjk[1]*uvkl[0];
	double hmm = (cross_ij_jk[0]*cross_jk_kl[0] + cross_ij_jk[1]*cross_jk_kl[1] + cross_ij_jk[2]*cross_jk_kl[2])/(sin(ang_ijk)*sin(ang_jkl));
	double tor;
	if(hmm < -1)
		tor = acos(-1)*180/acos(-1);
	else if(hmm > 1)
		tor = acos(1)*180/acos(-1);
	else
		tor = acos((cross_ij_jk[0]*cross_jk_kl[0] + cross_ij_jk[1]*cross_jk_kl[1] + cross_ij_jk[2]*cross_jk_kl[2])/(sin(ang_ijk)*sin(ang_jkl)))*180/acos(-1);
	double cross_x = cross_ij_jk[1] * cross_jk_kl[2] - cross_ij_jk[2] * cross_jk_kl[1];
	double cross_y = cross_ij_jk[2] * cross_jk_kl[0] - cross_ij_jk[0] * cross_jk_kl[2];
	double cross_z = cross_ij_jk[0] * cross_jk_kl[1] - cross_ij_jk[1] * cross_jk_kl[0];
	double norm = cross_x*cross_x + cross_y*cross_y + cross_z*cross_z;
	cross_x /= norm;
	cross_y /= norm;
	cross_z /= norm;
	double sign = 1.0;
	double dot = cross_x*uvjk[0] + cross_y*uvjk[1] + cross_z*uvjk[2];
	if(dot < 0.0) sign = -1.0;
	
	return tor*sign;
}

double Molecule::com(int cart)
{
	double cm = 0.0;
	double cmnum = 0.0;
	double cmden= 0.0;
	for(int i=0; i < natom; i++)
	{
		cmnum += geom[i][cart]*amasses[zvals[i]];
		cmden += amasses[zvals[i]];
	}
	cm = cmnum/cmden;
	return cm;			
}

void Molecule::translate(double x, double y, double z)
{
	for(int i=0; i < natom; i++)
	{
		geom[i][0] += x;
		geom[i][1] += y;
		geom[i][2] += z;
	}
}

void Molecule::moi()
{
	double mi;
	for(int k=0; k < natom; k++)
	{
		mi = amasses[zvals[k]];
		moit[0][0] += mi*(pow(geom[k][1],2.0) + pow(geom[k][2],2.0));
		moit[0][1] += mi*geom[k][0]*geom[k][1];
		moit[0][2] += mi*geom[k][0]*geom[k][2];
		moit[1][0] += mi*geom[k][1]*geom[k][0];
		moit[1][1] += mi*(pow(geom[k][0],2.0) + pow(geom[k][2],2.0));
		moit[1][2] += mi*geom[k][1]*geom[k][2];
		moit[2][0] += mi*geom[k][2]*geom[k][0];
		moit[2][1] += mi*geom[k][2]*geom[k][1];
		moit[2][2] += mi*(pow(geom[k][0],2.0) + pow(geom[k][1],2.0));
	}
}

double Molecule::rot_const(double inert)
{
	double h = 6.62607004E-34;
	double _pi = acos(-1);
	double rot = h/(8.0 * _pi * _pi * inert);
	rot /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
	return rot;
}

Molecule::Molecule(const char *filename, int q)
{
	//open file
	ifstream is(filename);
	assert(is.good());

	is >> natom;

	//allocate mem for zvals, geom and bonds, read file for zvals and geom values
	//I did loop fusion!
	zvals = new int[natom];
	geom = new double* [natom];
	bonds = new double* [natom];
	moit = new double* [3];
	for(int i=0; i < 3; i++)
		moit[i] = new double[3];
	for(int i=0; i < natom; i++)
	{
		geom[i] = new double[3];
		bonds[i] = new double[i];
		is >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];
	}

	is.close();	

	//Hessian part;
	ifstream iss("hessian.dat");
	assert(iss.good());

	int natom2;
	iss >> natom2;
	assert(natom2 == natom);

    H = new double* [natom2*3];
	for(int i=0; i < natom2*3; i++)
	{
		H[i] = new double[natom2*3];
		for(int j=0; j < natom2; j++)
			iss >> H[i][3*j] >> H[i][3*j+1] >> H[i][3*j+2];
	}

	iss.close();
}

Molecule::~Molecule()
{
	delete[] zvals;
	for(int i=0; i < natom; i++)
	{
		delete[] geom[i];
		delete[] bonds[i];
	}
	for(int i=0; i < 3; i++)
		delete[] moit[i];
	delete[] moit;
	delete[] geom;
	delete[] bonds;
	for(int i=0; i < natom*3; i++)
	{
		delete[] H[i];
	}
	delete[] H;
}
