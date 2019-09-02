// gpuswaxs.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#include <iostream>
#include <fstream>
#include <arrayfire.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <map>

using namespace af;

typedef unsigned char byte;

#include "aff.h"
#include "voxel.h"
#include "pdb.h"



int main()
{
    std::cout << "Hello World!\n";
	std::cout << __cplusplus << std::endl;

	for (int u = 0; u < 9; ++u) {
		std::cout << COEFS[H1n][u] << "  ";
	}
	std::cout << std::endl;

	//Voxel vox("por.binvox");
	//vox.swaxs(0.0, 0.002, 1.0, 0.5, 1.0);
	//std::cout << vox.getdim() << std::endl;

	PDB solute("solute.pdb");
	std::string* res = solute.get_residues();
	std::string* aname = solute.get_atomnames();
	/*
	for (int u = 0; u < solute.get_dim(); ++u) {
		std::cout << aname[u] << std::endl;
	}*/



	std::cout << "H id = " << solute.calculate_form_factor(AFMAP["H"], 1.0) << std::endl;
	

	return 0;
}








// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
