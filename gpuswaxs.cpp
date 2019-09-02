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
#include <set>


using namespace af;
typedef unsigned char byte;

#include "cxxopts.hpp"

#include "aff.h"
#include "voxel.h"
#include "pdb.h"



// The parsing function
cxxopts::ParseResult parse(int argc, char* argv[]){
	
	try {

		cxxopts::Options options("GPUSWAXS", "GPU implementation of SWAXS profile calculation. \n");
		options.positional_help("[ARGS...]").show_positional_help();

		options.add_options()
			("B,binvox", "The shape file of dummy voxels: .binvox", cxxopts::value<std::string>()->default_value(""))
			("P,pdb", "The single PDB mode for atomic coordinates: .pdb", cxxopts::value<std::string>()->default_value(""))
			("T,solute", "The solute.pdb file", cxxopts::value<std::string>()->default_value(""))
			("V,solvent", "The solvent.pdb file", cxxopts::value<std::string>()->default_value(""))
			("d,density", "The electron density averaged on a voxel, in e/A^3", cxxopts::value<std::string>()->default_value("0.50"))
			("v,voxdim", "The voxel size in A", cxxopts::value<std::string>()->default_value("2.0"))
			("J", "The number of orientations", cxxopts::value<std::string>()->default_value("1800"))
			("qmin", "The minimum q in 1/A", cxxopts::value<float>())
			("qspacing", "The spacing in q axis", cxxopts::value<float>())
			("qmax", "The maximum q reached", cxxopts::value<float>())
			("h,help", "Print this help message and exit");

		
		options.parse_positional({ "qmin", "qspacing", "qmax" });

		auto result = options.parse(argc, argv);


		if (result.count("help")) {
			std::cout << options.help({ "" }) << std::endl;
			exit(0);
		}

		// make sure that positional arguments are good. 
		if (result.count("qmin") == 0) {
			std::cout << "GPUSWAXS: The positional argument [qmin] is missing, try --help." << std::endl;
			exit(0);
		}

		if (result.count("qspacing") == 0) {
			std::cout << "GPUSWAXS: The positional argument [qspacing] is missing, try --help." << std::endl;
			exit(0);
		}

		if (result.count("qmax") == 0) {
			std::cout << "GPUSWAXS: The positional argument [qmax] is missing, try --help." << std::endl;
			exit(0); 
		}

		return result;
	}
	catch (const cxxopts::OptionException& e) {
		std::cout << "GPUSWAXS: Error parsing options: " << e.what() << std::endl;
		exit(1);
	}
}





int main(int argc, char* argv[])
{

	// command line header
	{
		std::cout << std::endl << std::endl;
		std::cout << std::string(80, '*') << std::endl;
		std::cout << "**********     gpuswaxs.exe: version 1.0.0" << std::string(80 - 27 - 25, ' ') << "**********" << std::endl;
		std::cout << "**********     A GPU-accelerated computation of the small-wide-angle" << std::string(80 - 53 - 25, ' ') << "**********" << std::endl;
		std::cout << "**********     X-ray scattering (SWAXS) profile (.dat) from .binvox " << std::string(80 - 53 - 25, ' ') << "**********" << std::endl;
		std::cout << "**********     or .pdb file(s). Require: ArrayFire Library" << std::string(80 - 43 - 25, ' ') << "**********" << std::endl;
		std::cout << "**********" << std::string(80 - 20, ' ') << "**********" << std::endl;

		std::cout << "**********     AUTHOR: Yen-Lin Chen @ Cornell" << std::string(80 - 30 - 25, ' ') << "**********" << std::endl;
		std::cout << "**********     DATE:   Sept. 2 2019" << std::string(80 - 20 - 25, ' ') << "**********" << std::endl;
		std::cout << std::string(80, '*') << std::endl << std::endl;
	}



	auto parsed = parse(argc, argv);


	// parse all the options
	std::string binvox = parsed["binvox"].as<std::string>();
	std::string pdb = parsed["pdb"].as<std::string>();
	std::string solute = parsed["solute"].as<std::string>();
	std::string solvent = parsed["solvent"].as<std::string>();
	float density = std::stof(parsed["density"].as<std::string>());
	float voxdim = std::stof(parsed["voxdim"].as<std::string>());
	int J = std::stoi(parsed["J"].as<std::string>());
	float qmin = parsed["qmin"].as<float>();
	float qspacing = parsed["qspacing"].as<float>();
	float qmax = parsed["qmax"].as<float>();


	// Check the values of the arguments
	if (qspacing <= 0.0) {
		std::cout << "GPUSWAXS: qspacing should be positive, but qspacing = " << qspacing << " detected ..." << std::endl;
		exit(1); 
	}

	if (qmin < 0.0) {
		std::cout << "GPUSWAXS: qmin should be non-negative, but qspacing = " << qspacing << " detected ..." << std::endl;
		exit(1);
	}

	if (qmin >= qmax) {
		std::cout << "GPUSWAXS: qmin should be less than qmax, but the specified values: qmin = " << qmin << " > qmax = " << qmax << " ..." << std::endl;
		exit(1);
	}


	// Do the job based on binvox, pdb, solute, solvent
	if ((binvox != "") && (pdb == "") && (solute == "") && (solvent == "")) {
		try {

			if (density < 0) {
				std::cout << "GPUSWAXS: Negative density = " << density << " detected ..." << std::endl;
				exit(1);
			}

			if (voxdim < 0) {
				std::cout << "GPUSWAXS: Negative voxdim = " << voxdim << " detected ..." << std::endl;
				exit(1);
			}

			std::cout << "GPUSWAXS: Operating in binvox shape mode ..." << std::endl;
			std::cout << "GPUSWAXS: binvox file = " << binvox << " with J = " << J << " and q = " << qmin << ":" << qspacing << ":" << qmax << " ";
			std::cout << "density = " << density << " voxdim = " << voxdim << std::endl; 
			Voxel vox(binvox);
			vox.swaxs(qmin, qspacing, qmax, density, pow(voxdim, 3), J);
		}
		catch (...) {
			std::cout << "GPUSWAXS: Error encountered. Please check --help or file issue. " << std::endl;
		}
	}

	else if ((binvox == "") && (pdb != "") && (solute == "") && (solvent == "")) {
		try {
			std::cout << "GPUSWAXS: Operating in single pdb mode ..." << std::endl;
			std::cout << "GPUSWAXS: pdb file = " << pdb << " with J = " << J << " and q = " << qmin << ":" << qspacing << ":" << qmax << " " << std::endl;
			PDB pdb1(pdb);
			pdb1.swaxs(qmin, qspacing, qmax, J);
		}
		catch (...) {
			std::cout << "GPUSWAXS: Error encountered. Please check --help or file issue. " << std::endl;
		}
	}

	else if ((binvox == "") && (pdb == "") && (solute != "") && (solvent != "")) {
		try {
			std::cout << "GPUSWAXS: Operating in pair pdb mode ..." << std::endl;
			std::cout << "GPUSWAXS: solute file = " << solute << " solvent file = " << solvent << " with J = " << J << " and q = " << qmin << ":" << qspacing << ":" << qmax << " " << std::endl;
			PDB pdbsolute(solute);
			PDB pdbsolvent(solvent); 
			pdbsolute.swaxs(pdbsolvent, qmin, qspacing, qmax, J);
		}
		catch (...) {
			std::cout << "GPUSWAXS: Error encountered. Please check --help or file issue. " << std::endl;
		}
	}
	else {
		std::cout << "GPUSWAXS: There might be conflicting options ..." << std::endl;
		std::cout << "GPUSWAXS: Check the following arguments to make sure of that ... " << std::endl;
		std::cout << "          1. You specified either [binvox], [pdb] or [solute, solvent]. " << std::endl;
		std::cout << "          2. [solute, solvent] should be specified or empty at the same time. " << std::endl << std::endl;
		std::cout << "GPUSWAXS: binvox = [" << binvox << "] pdb = [" << pdb << "] solute = [" << solute << "] solvent = [" << solvent << "]" << std::endl;
	}

	
	return EXIT_SUCCESS;
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
