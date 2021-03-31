#include <fstream>
#include <cmath>
#include <cassert>
#include "hf.h"
#include "masses.h"
#include <iostream>
#include <sstream>
#include <string>


using namespace std;

void hf::example()
{
}

//Creation
hf::hf(const char *mol, const char *basis, int q)
{
	string a = mol;
	string b = basis;
	string s = "input/" + a + "/" + b + "/";
	ifstream is_enuc(s + "enuc.dat");
	assert(is_enuc.good());

	is_enuc >> enuc;
	is_enuc.close();
}

//Annhiliation
hf::~hf()
{

}
