#pragma once

class Voxel {

private:

	// private members
	int nr_voxels;
	array coordata;
	std::string fn;

public:

	// constructors
	Voxel(int n); 
	Voxel(std::string filename); 

	// public methods
	array getdata() const { return coordata; };
	int getdim() const { return nr_voxels; };
	std::string get_filename() const { return fn; };
	void swaxs(float qmin, float qspacing, float qmax, float d, float voxvol, int J);
};



// null constructors
Voxel::Voxel(int n) {
	nr_voxels = n;
	coordata = array(n, 3, s32);
	fn = "empty";
}



// main constructor
Voxel::Voxel(std::string filename) {

	fn = filename;

	std::ifstream* input = new std::ifstream(filename.c_str(), std::ios::in | std::ios::binary);

	if (!input->good()) {
		std::cout << "Voxel Error: Cannot open [ " << fn << " ]. " << std::endl;
		exit(1);
	}

	// read binvox header
	std::string line;
	*input >> line;
	if (line.compare("#binvox") != 0) {
		std::cout << "Voxel Error: First line reads [" << line << "] instead of [#binvox]" << std::endl;
		delete input;
		Voxel(0);
	}
	else {
		std::cout << "Voxel INFO: [#binvox]" << std::endl;
	}

	int version;
	*input >> version;
	std::cout << "Voxel INFO: Reading binvox version: " << version << std::endl;

	// binvox dimensions
	int depth = -1, height = 0, width = 0;
	int done = 0;
	while (input->good() && !done) {
		*input >> line;
		if (line.compare("data") == 0) done = 1;
		else if (line.compare("dim") == 0) {
			*input >> depth >> height >> width;
		}
		else {
			// Ignore the translate and scale
			std::cout << "Voxel INFO: Keyword [" << line << "], ignored." << std::endl;
			char c;
			do {
				c = input->get();
			} while (input->good() && (c != '\n'));
		}
	}

	if (!done) {
		std::cout << "Voxel Error: Cannot readin the headers." << std::endl;
		Voxel(0);
	}

	if (depth == -1) {
		std::cout << "Voxel Error: Missing dimensions in the headers." << std::endl;
		Voxel(0);
	}

	// Prepare to read in the binary data
	int size = width * height * depth;
	byte* voxels = new byte[size];
	if (!voxels) {
		std::cout << "Voxel Error: Cannot allocating memory." << std::endl;
		Voxel(0);
	}

	// read binvox data
	byte value, count;
	int index = 0, end_index = 0;
	nr_voxels = 0; 

	input->unsetf(std::ios::skipws);
	*input >> value;

	while ((end_index < size) && input->good()) {
		*input >> value >> count;

		if (input->good()) {
			end_index = index + count;
			if (end_index > size) {
				std::cout << "Voxel Error: Unexpected voxel reading error." << std::endl;
				coordata = array(0, 3, s32);
				nr_voxels = 0;
			}
			for (int i = index; i < end_index; i++) voxels[i] = value;

			if (value) nr_voxels += count;
			index = end_index;
		}
	}

	input->close();
	std::cout << "Voxel INFO: Success, read " << nr_voxels << " voxels" << std::endl;

	// Take care of the zero-voxel, will cause error later
	if (nr_voxels == 0) exit(0);

	// Use a long 1D array to store the x, y, z coordinates to be mapped into GPU later
	// structure: [x1, x2 ... xn, y1, y2 ... yn, z1, z2 ... zn]
	int* a = new int[3 * nr_voxels];

	int m = 0;

	for (int i = 0; i < size; i++) {
		if (int(voxels[i]) == 0) continue;
		else {
			a[m] = i / (width * height);
			int ii = i % (width * height);
			a[m + nr_voxels] = ii / height;
			a[m + 2 * nr_voxels] = ii % height;
			m++;
		}
	}

	// copy a to GPU
	array tmp(nr_voxels, 3, a);
	coordata = 2.0 * tmp;

	delete[] a;

	std::cout << "Voxel INFO: " << filename << " read-in successfully." << std::endl;
}




void Voxel::swaxs(float qmin, float qspacing, float qmax, float d, float voxvol=8.0, int J=2500) {


	std::string filename = fn;
	int nq, i, j, m = getdim();
	nq = static_cast<int>(ceil((qmax - qmin) / qspacing));

	// Take care of the matmul-error from empty matrix 
	if (m == 0) {
		std::cout << "Voxel INFO: No voxel found in [ " << fn << " ], program terminating .." << std::endl;
		exit(0);
	}
	float* q = new float[nq];
	float* intensity = new float[nq];
	for (i = 0; i < nq; i++) q[i] = qmin + qspacing * i;

	// solid angles
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

	array qr(m, J, f32);
	array sys1(2, J, f32), amplitude(1, J, f32);
	
	std::cout << "Voxel INFO: Start to calculate swaxs curves with J = " << J << " and q = " << qmin << ":" << qspacing << ":" << qmax << " density = " << d * voxvol << "e/A^3" << std::endl;
	timer::start();

	// matrix operations using GPU for every q[i]
	for (i = 0; i < nq; i++) {
		qr = (d * voxvol) * q[i] * matmul(coordata, qmat);
		sys1.row(0) = sum(cos(qr));
		sys1.row(1) = sum(-sin(qr)); 

		amplitude = sum(pow(sys1, 2));
		intensity[i] = mean(amplitude, 1)(0).scalar<float>();

		// cool loading bar
		std::cout << "\r" << (100 * i / (nq - 1)) << "% completed: ";
		std::cout << std::string((50 * i / nq), '|');
		std::cout.flush();
	}

	printf("\nVoxel INFO: Elapsed time: %g seconds\n", timer::stop());

	// print out the data file
	filename.replace(filename.end() - 6, filename.end(), "dat");
	std::ofstream* out = new std::ofstream(filename);
	if (!out->good()) {
		std::cout << "Voxel Error: Cannot open [ " << filename << " ]. " << std::endl;
		exit(1);
	}

	std::cout << "Voxel INFO: Writing q, and intensity to the file: " << filename << " ." << std::endl;
	(*out).precision(17);
	for (int i = 0; i < nq; i++) {
		*out << q[i] << "\t" << intensity[i] << std::endl;
	}
	out->close();
	std::cout << "Voxel INFO: Writing completed successfully." << std::endl;

	delete[] q;
	delete[] intensity;
}