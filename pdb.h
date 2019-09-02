#pragma once

/*------------------------class PDB--------------------------------------*/
class PDB {

private:
	std::string fn; // filename
	float* coordinates; // atomic coordinates
	int n; // number of atoms in the pdb
	std::string* atomnames;
	std::string* residues;

public:
	// Constructors
	PDB(int m);
	PDB(std::string filename);

	// Functions
	std::string get_filename() const { return fn; };
	float* get_data() const { return coordinates; };
	int get_dim() const { return n; };
	std::string* get_atomnames() const { return atomnames; };
	std::string* get_residues() const { return residues; };
	int atom_match(std::string str) const; // process the atomnames 
	float calculate_form_factor(int id, float q) const; // Calculate a vector of form factor at q 
	void swaxs(PDB solvent, float qmin, float qspacing, float qmax, int J); // Calculate swaxs curve using solute and solvent pdbs
	void swaxs(float qmin, float qspacing, float qmax, int J); // Calculate swaxs curve using one pdb (in vacuo or in vitro)

};




// null constructor
PDB::PDB(int m) {
	n = m;
	fn = "empty";
	coordinates = new float[3 * n];
	atomnames = new std::string[n];
	residues = new std::string[n];
}




// Main Constructor
PDB::PDB(std::string filename) {
	// Read in the pdb files
	fn = filename;

	std::ifstream input(fn.c_str(), std::ios::in);

	if (input.fail()) {
		std::cout << "Error: Cannot open [ " << fn << " ]. " << std::endl;
		exit(1);
	}

	std::string line;
	int natoms = 0;

	// First Scan and check the number of atoms in order to use static arrays
	while (!input.eof()) {
		std::getline(input, line);
		if (line.find("ATOM") == 0) natoms++;
	}
	input.close();
	std::cout << "INFO: Detected " << natoms << " atoms in the file: " << fn << "." << std::endl;
	// if (natoms > 50000) std::cout << "WARN: There are too many atoms in the pdb file.\n      Consider using parallel Julia. " << std::endl;
	n = natoms;


	// read again and store relevant information using static array
	std::cout << "INFO: Processing all atomic coordinates ... " << std::endl;
	input = std::ifstream(fn.c_str(), std::ios::in);
	coordinates = new float[3 * n];
	atomnames = new std::string[n];
	residues = new std::string[n];
	natoms = 0;
	std::string* elements = new std::string[13];
	int atom_col = 10;


	while (!input.eof()) {
		std::getline(input, line);
		if (line.find("ATOM") == 0) {
			int col = 0;
			std::istringstream iss(line);

			// Parse pdb columns and locate the last meaningful column
			for (int i = 0; i < 13; i++) {
				iss >> elements[i];
				if (elements[i].length() != 0) atom_col = i;
			}


			// residues are always the 4th column
			residues[natoms] = elements[3];

			// coordinates and atomnames should be traced from the last column
			atomnames[natoms] = elements[atom_col];
			coordinates[natoms] = std::stof(elements[atom_col - 5].c_str());
			coordinates[natoms + n] = std::stof(elements[atom_col - 4].c_str());
			coordinates[natoms + 2 * n] = std::stof(elements[atom_col - 3].c_str());

			natoms++;

		}
	}
	delete[] elements;

	if (natoms != n) std::cout << "WARN: Atom number mismatch." << std::endl;

	// Convert atom names: for SOL-O, SOL-H, my convention to take care of the electron withdrawing effects
	for (int i = 0; i < n; i++) {
		if ((residues[i] == "WAT") || (residues[i] == "SOL") || (residues[i] == "HOH")) {
			atomnames[i].replace(0, 0, "SOL");
		}
	}

	std::cout << "INFO: Finished processing all atomic coordinates. " << std::endl;
}




// Calculate unique form factors for atom id at a specific q value
float PDB::calculate_form_factor(int id, float q) const {
	float sum = COEFS[id][8];
	for (int j = 0; j < 4; j++) {
		sum = sum + COEFS[id][j] * exp(-COEFS[id][j + 4] * pow((0.25 * q / Pi), 2));
	}
	return sum;
}



// 





/*
// Calculate swaxs curve in vacuo or in vitro.. (the shape of the solvent contributes if solvent exists... )
// NOTE: This is probably more useful for placevent-ed conformations without bulk solvent..
void PDB::swaxs() {

	std::string filename = fn;

	// Initialize q and intensity
	int nq = 701, i, j, J = 1800, m = n;
	float* q = new float[nq];
	float* intensity = new float[nq];
	for (i = 0; i < nq; i++) q[i] = 0.002 * i;


	// Transfer PDB data to GPU
	array mat(m, 3, (float*)coordinates);
	array qr(m, J, f32);
	array re(1, J, f32), im(1, J, f32);


	// Solid angle average
	float* xx = new float[J];
	array theta(1, J, f32), phi(1, J, f32);
	for (j = 0; j < J; j++)
		xx[j] = (2.0 * (j + 1.0) - 1.0 - J) / J;

	array x(1, J, xx);
	delete[] xx;
	theta = acos(x);
	phi = sqrt(Pi * J) * asin(x);

	array qmat = randu(3, J, f32);
	qmat.row(0) = sin(theta) * cos(phi);
	qmat.row(1) = sin(theta) * sin(phi);
	qmat.row(2) = cos(theta);

	// allocate memory for the form factor vector
	float* atomff = new float[m];

	std::cout << "INFO: Start to calculate swaxs curves..." << std::endl;
	timer::start();

	// Matrix operations...
	for (int k = 0; k < nq; k++) {

		// Get a list of form factors 
		double* ff = calculate_form_factor(q[k]);

		for (int i = 0; i < m; i++) {
			atomff[i] = ff[atom_match(atomnames[i])];
		}

		array af(1, m, (float*)atomff);

		qr = q[k] * matmul(mat, qmat);
		re.row(0) = sum(matmul(af, cos(qr)), 0);
		im.row(0) = sum(matmul(af, sin(qr)), 0);
		im = pow(im, 2) + pow(re, 2);
		intensity[k] = mean(im)(0).scalar<float>();

		// Cool Loading bar
		std::cout << "\r" << (100 * k / (nq - 1)) << "% completed: ";
		std::cout << std::string((50 * k / nq), '|');
		std::cout.flush();
	}
	printf("\nINFO: Elapsed time: %g seconds\n", timer::stop());

	// The output scattering curve is not normalized, the unit is in e^2
	filename.replace(filename.end() - 3, filename.end(), "dat");
	std::ofstream* out = new std::ofstream(filename);
	if (!out->good()) {
		std::cout << "Error: Cannot open [ " << filename << " ]. " << std::endl;
		exit(1);
	}

	std::cout << "INFO: Writing q, and intensity profile to the file: " << filename << " ." << std::endl;
	for (int k = 0; k < nq; k++) {
		*out << q[k] << "\t" << intensity[k] << std::endl;
	}
	out->close();
	std::cout << "INFO: Writing completed successfully." << std::endl;

	delete[] atomff;
	delete[] q;
	delete[] intensity;
}

// Solute and solvent subtraction using Park. et al
void PDB::swaxs(PDB solvent)
{

	std::string filename = fn;

	// Initialize q and intensity
	int nq = 701, i, j, J = 1800, m1 = n, m2 = solvent.get_dim();
	float* q = new float[nq];
	float* intensity = new float[nq];
	for (i = 0; i < nq; i++) q[i] = 0.002 * i;

	// Transfer PDB data to GPU
	//for (int i = 0; i < 3 * n; i++) std::cout << coordinates[i] << " ";
	array mat1(m1, 3, (float*)coordinates);
	array mat2(m2, 3, (float*)solvent.get_data());
	array qr1(m1, J, f32);
	array qr2(m2, J, f32);
	array re(1, J, f32), im(1, J, f32);

	// Solid angle average
	float* xx = new float[J];
	array theta(1, J, f32), phi(1, J, f32);
	for (j = 0; j < J; j++)
		xx[j] = (2.0 * (j + 1.0) - 1.0 - J) / J;

	array x(1, J, xx);
	delete[] xx;
	theta = acos(x);
	phi = sqrt(Pi * J) * asin(x);

	array qmat = randu(3, J, f32);
	qmat.row(0) = sin(theta) * cos(phi);
	qmat.row(1) = sin(theta) * sin(phi);
	qmat.row(2) = cos(theta);

	float* atomff1 = new float[m1];
	float* atomff2 = new float[m2];

	std::cout << "INFO: Start to calculate swaxs curves..." << std::endl;
	timer::start();

	// matrix operations
	for (int k = 0; k < nq; k++) {

		// Get a list of form factors 
		double* ff1 = calculate_form_factor(q[k]);
		double* ff2 = solvent.calculate_form_factor(q[k]);
		for (int i = 0; i < m1; i++) atomff1[i] = ff1[atom_match(atomnames[i])];
		for (int i = 0; i < m2; i++) atomff2[i] = ff2[atom_match(atomnames[i])];

		array af1(1, m1, (float*)atomff1), af2(1, m2, (float*)atomff2);

		qr1 = q[k] * matmul(mat1, qmat);
		qr2 = q[k] * matmul(mat2, qmat);
		re = matmul(af1, cos(qr1)) - matmul(af2, cos(qr2));
		// minus sign here matters due to different atomic form factors 
		im = matmul(af1, -sin(qr1)) - matmul(af2, -sin(qr2));
		im = pow(im, 2) + pow(re, 2);

		intensity[k] = mean(im)(0).scalar<float>();

		// Cool Loading bar
		std::cout << "\r" << (100 * k / (nq - 1)) << "% completed: ";
		std::cout << std::string((50 * k / nq), '|');
		std::cout.flush();
	}
	printf("\nINFO: Elapsed time: %g seconds\n", timer::stop());

	// save file, again, the filename is changed too. 
	filename.replace(filename.end() - 4, filename.end(), "_buffer_subtracted.dat");
	std::ofstream* out = new std::ofstream(filename);
	if (!out->good()) {
		std::cout << "Error: Cannot open [ " << filename << " ]. " << std::endl;
		exit(1);
	}

	std::cout << "INFO: Writing q, and intensity to the file: " << filename << " ." << std::endl;
	for (int k = 0; k < nq; k++) {
		*out << q[k] << "\t" << intensity[k] << std::endl;
	}
	out->close();
	std::cout << "INFO: Writing completed successfully." << std::endl;

	delete[] atomff1;
	delete[] atomff2;
	delete[] q;
	delete[] intensity;

}
*/

