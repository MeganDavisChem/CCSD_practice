#include <string>

using namespace std;

class Molecule
{
	public:
		int natom;
		int charge;
		int *zvals;
		double **geom;
		double **bonds;
		double **moit;
	
		void print_geom();
		void translate(double x, double y, double z);
		double bond_lengths(int atom1, int atom2);
		double bond_angles(int atom1, int atom2, int atom3);
		double unit_vectors(int atom1, int atom2, int cart);
		double oop_angles(int atom1, int atom2, int atom3, int atom4);
		double torsion(int atom1, int atom2, int atom3, int atom4);
		double com(int cart);
		double rot_const(double inert);
		void moi();
	
		Molecule(const char *filename, int q);
		~Molecule();
};

