#pragma once

/*------------------------class PDB--------------------------------------*/
class PDB {

private:
	std::string fn; // filename
	float* coordinates; // atomic coordinates
	int n; // number of atoms in the pdb
	std::string* atomnames;
	std::string* residues;
	float calculate_form_factor(int id, float q) const; // Calculate a vector of form factor at q 


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
	std::map<std::string, float> uniquemap(float q);
	void swaxs(float qmin, float qspacing, float qmax, int J); // Calculate swaxs curve using one pdb (in vacuo or in vitro)
	void swaxs(PDB solvent, float qmin, float qspacing, float qmax, int J); // Calculate swaxs curve using solute and solvent pdbs


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
		std::cout << "PDB Error: Cannot open [ " << fn << " ]. " << std::endl;
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
	std::cout << "PDB INFO: Detected " << natoms << " atoms in the file: " << fn << "." << std::endl;
	// if (natoms > 50000) std::cout << "WARN: There are too many atoms in the pdb file.\n      Consider using parallel Julia. " << std::endl;
	n = natoms;


	// read again and store relevant information using static array
	std::cout << "PDB INFO: Processing all atomic coordinates ... " << std::endl;
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

	if (natoms != n) std::cout << "PDB WARN: Atom number mismatch." << std::endl;

	// Convert atom names: for SOL-O, SOL-H, my convention to take care of the electron withdrawing effects
	for (int i = 0; i < n; i++) {
		if ((residues[i] == "WAT") || (residues[i] == "SOL") || (residues[i] == "HOH")) {
			atomnames[i].replace(0, 0, "SOL");
		}
	}

	std::cout << "PDB INFO: Finished processing all atomic coordinates. " << std::endl;
}




// Calculate unique form factors for atom id at a specific q value
float PDB::calculate_form_factor(int id, float q) const {
	float sum = COEFS[id][8];
	for (int j = 0; j < 4; j++) {
		sum = sum + COEFS[id][j] * exp(-COEFS[id][j + 4] * pow((0.25 * q / Pi), 2));
	}
	return sum;
}




// calculate a unique map at a specifit q
std::map<std::string, float> PDB::uniquemap(float q) {

	std::set<std::string> unique(atomnames, atomnames + n); 
	std::set<std::string>::iterator p_unique = unique.begin();

	int nunique = unique.size();
	std::map<std::string, float> rlt; 
	std::string tmp; 
	for (int i = 0; i < nunique; ++i) {
		tmp = *p_unique;
		rlt.insert(std::pair<std::string, float>(tmp, calculate_form_factor(AFMAP[tmp], q)));
		p_unique++;
	}


	// Deal with the SOLO and SOLH
	if (rlt.count("SOLH") == 1) rlt.at("SOLH") = rlt["SOLH"] * (1.0 - 0.48 * exp(-0.5 * pow(q / 2.2, 2)));
	if (rlt.count("SOLO") == 1) rlt.at("SOLO") = rlt["SOLO"] * (1.0 + 0.12 * exp(-0.5 * pow(q / 2.2, 2)));

	return rlt; 
	
}






// Calculate swaxs curve in vacuo or in vitro.. (the shape of the solvent contributes if solvent exists... )
// NOTE: This is probably more useful for placevent-ed conformations without bulk solvent..
void PDB::swaxs(float qmin, float qspacing, float qmax, int J) {


	std::string filename = fn;

	// Initialize q and intensity
	int nq, i, j, m = get_dim();
	nq = static_cast<int>(ceil((qmax - qmin) / qspacing));
	
	float* q = new float[nq];
	float* intensity = new float[nq];
	for (i = 0; i < nq; i++) q[i] = qmin + qspacing * i;


	// Transfer PDB data to GPU
	array mat(m, 3, coordinates);


	// Solid angle average
	float* xx = new float[J];
	array theta(1, J, f32), phi(1, J, f32);
	for (j = 0; j < J; j++) xx[j] = (2.0 * (j + 1.0) - 1.0 - J) / J;

	array x(1, J, xx);
	delete[] xx;
	theta = acos(x);
	phi = sqrt(Pi * J) * asin(x);

	array qmat(3, J, f32);
	qmat.row(0) = sin(theta) * cos(phi);
	qmat.row(1) = sin(theta) * sin(phi);
	qmat.row(2) = cos(theta);

	array qr = matmul(mat, qmat);
	array sys1(2, J, f32), amplitude(1, J, f32);



	// allocate memory for the form factor vector
	float* atomff = new float[m];

	std::cout << "PDB INFO: Start to calculate swaxs curves..." << std::endl;
	timer::start();

	// Matrix operations...
	for (int k = 0; k < nq; k++) {

		// Get a list of form factors 
		std::map<std::string, float> atommap = uniquemap(q[k]);
		for (i = 0; i < m; i++) atomff[i] = atommap[atomnames[i]];
		array aff(1, m, atomff);

		sys1.row(0) = matmul(aff, cos(q[k] * qr));
		sys1.row(1) = matmul(aff, -sin(q[k] * qr));
		amplitude = sum(pow(sys1, 2));
		intensity[k] = mean(amplitude)(0).scalar<float>();

		// Cool Loading bar
		std::cout << "\r" << (100 * k / (nq - 1)) << "% completed: ";
		std::cout << std::string((50 * k / nq), '|');
		std::cout.flush();
	}
	printf("\nPDB INFO: Elapsed time: %g seconds\n", timer::stop());

	

	filename.replace(filename.end() - 3, filename.end(), "dat");
	std::ofstream* out = new std::ofstream(filename);
	if (!out->good()) {
		std::cout << "PDB Error: Cannot open [ " << filename << " ]. " << std::endl;
		exit(1);
	}

	std::cout << "PDB INFO: Writing q, and intensity profile to the file: " << filename << " ." << std::endl;
	for (int k = 0; k < nq; k++) {
		*out << q[k] << "\t" << intensity[k] << std::endl;
	}
	out->close();
	std::cout << "PDB INFO: Writing completed successfully." << std::endl;

	delete[] atomff;
	delete[] q;
	delete[] intensity;
}






// Solute and solvent subtraction using Park. et al
void PDB::swaxs(PDB solvent, float qmin, float qspacing, float qmax, int J)
{

	std::string filename = fn;
	std::string* atomnames2 = solvent.get_atomnames();

	// Initialize q and intensity
	int nq, i, j, m1 = get_dim(), m2 = solvent.get_dim();
	nq = static_cast<int>(ceil((qmax - qmin) / qspacing));

	float* q = new float[nq];
	float* intensity = new float[nq];
	for (i = 0; i < nq; i++) q[i] = qmin + qspacing * i;

	// Transfer PDB data to GPU
	//for (int i = 0; i < 3 * n; i++) std::cout << coordinates[i] << " ";
	array mat1(m1, 3, coordinates);
	array mat2(m2, 3, solvent.get_data());


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


	array qr1 = matmul(mat1, qmat);
	array qr2 = matmul(mat2, qmat);
	array sys1(2, J, f32), sys2(2, J, f32), amplitude(1, J, f32);

	float* atomff1 = new float[m1];
	float* atomff2 = new float[m2];

	std::cout << "PDB INFO: Start to calculate swaxs curves..." << std::endl;
	timer::start();

	// matrix operations
	for (int k = 0; k < nq; k++) {

		// Get a list of form factors 
		std::map<std::string, float> atommap1 = uniquemap(q[k]);
		std::map<std::string, float> atommap2 = solvent.uniquemap(q[k]);
		for (i = 0; i < m1; ++i) atomff1[i] = atommap1[atomnames[i]];
		for (i = 0; i < m2; ++i) atomff2[i] = atommap2[atomnames2[i]];
		array aff1(1, m1, atomff1);
		array aff2(1, m2, atomff2);

		sys1.row(0) = matmul(aff1, cos(q[k] * qr1));
		sys1.row(1) = matmul(aff1, -sin(q[k] * qr1));
		sys2.row(0) = matmul(aff2, cos(q[k] * qr2));
		sys2.row(1) = matmul(aff2, -sin(q[k] * qr2));

		amplitude = sum(pow(sys1 - sys2, 2));
		intensity[k] = mean(amplitude)(0).scalar<float>();

		// Cool Loading bar
		std::cout << "\r" << (100 * k / (nq - 1)) << "% completed: ";
		std::cout << std::string((50 * k / nq), '|');
		std::cout.flush();
	}
	printf("\nPDB INFO: Elapsed time: %g seconds\n", timer::stop());

	// save file, again, the filename is changed too. 
	filename.replace(filename.end() - 4, filename.end(), "_buffer_subtracted.dat");
	std::ofstream* out = new std::ofstream(filename);
	if (!out->good()) {
		std::cout << "PDB Error: Cannot open [ " << filename << " ]. " << std::endl;
		exit(1);
	}

	std::cout << "PDB INFO: Writing q, and intensity to the file: " << filename << " ." << std::endl;
	for (int k = 0; k < nq; k++) {
		*out << q[k] << "\t" << intensity[k] << std::endl;
	}
	out->close();
	std::cout << "PDB INFO: Writing completed successfully." << std::endl;

	delete[] atomff1;
	delete[] atomff2;
	delete[] q;
	delete[] intensity;

}


