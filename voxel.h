#pragma once

/*------------------------class Voxel--------------------------------------*/
class Voxel {
private:
	int nr_voxels; // number of voxels
	array coordata; // all the coordinates (x,y,z) of voxels stored in GPU array
	std::string fn; // filename
public:
	Voxel(int n); // null constructor
	Voxel(std::string filename); // The meaningful constructor from a .binvox filename
	array getdata() const { return coordata; };
	int getdim() const { return nr_voxels; };
	std::string get_filename() const { return fn; };
	void swaxs(); // Calculate SWAXS curevs and outputs the filename.dat
};

Voxel::Voxel(int n) {
	nr_voxels = n;
	coordata = array(n, 3, s32);
	fn = "empty";
}

// Construct Voxel by reading the .binvox file
Voxel::Voxel(std::string filename) {

	fn = filename;

	// input file stream pointer, use "->function()" as (input*).function()
	std::ifstream* input = new std::ifstream(filename.c_str(), std::ios::in | std::ios::binary);

	if (!input->good()) {
		std::cout << "Error: Cannot open [ " << fn << " ]. " << std::endl;
		exit(1);
	}

	// read binvox header
	std::string line;
	*input >> line;
	if (line.compare("#binvox") != 0) {
		std::cout << "Error: First line reads [" << line << "] instead of [#binvox]" << std::endl;
		delete input;
		Voxel(0);
	}
	else {
		std::cout << "INFO: [#binvox]" << std::endl;
	}

	int version;
	*input >> version;
	std::cout << "INFO: Reading binvox version: " << version << std::endl;

	// binvox dimensions
	int depth, height, width;
	depth = -1;
	int done = 0;
	while (input->good() && !done) {
		*input >> line;
		if (line.compare("data") == 0) done = 1;
		else if (line.compare("dim") == 0) {
			*input >> depth >> height >> width;
		}
		else {
			// Ignore the translate and scale
			std::cout << "INFO: Keyword [" << line << "], ignored." << std::endl;
			char c;
			do {
				c = input->get();
			} while (input->good() && (c != '\n'));
		}
	}

	if (!done) {
		std::cout << "Error: Cannot readin the headers." << std::endl;
		Voxel(0);
	}

	if (depth == -1) {
		std::cout << "Error: Missing dimensions in the headers." << std::endl;
		Voxel(0);
	}

	// Prepare to read in the binary data
	int size = width * height * depth;
	byte* voxels = new byte[size];
	if (!voxels) {
		std::cout << "Error: Cannot allocating memory." << std::endl;
		Voxel(0);
	}

	// read binvox data
	byte value, count;
	int index = 0, end_index = 0;
	nr_voxels = 0; // It is built in instance: this.nr_voxels

	input->unsetf(std::ios::skipws);
	*input >> value;

	while ((end_index < size) && input->good()) {
		*input >> value >> count;

		if (input->good()) {
			end_index = index + count;
			if (end_index > size) {
				std::cout << "Error: Unexpected voxel reading error." << std::endl;
				coordata = array(0, 3, s32);
				nr_voxels = 0;
			}
			for (int i = index; i < end_index; i++) voxels[i] = value;

			if (value) nr_voxels += count;
			index = end_index;
		}
	}

	input->close();
	std::cout << "INFO: Success, read " << nr_voxels << " voxels" << std::endl;

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

	// Copy a to GPU
	array tmp(nr_voxels, 3, a);
	coordata = 2.0 * tmp;

	delete[] a;

	std::cout << "INFO: " << filename << " read-in successfully." << std::endl;
}

void Voxel::swaxs() {

	// NOTE: use the type of f32 instead of f64 and the performance is improved by 3x

	std::string filename = fn;
	int nq = 701, i, j, J = 2500, m = getdim();

	// Take care of the matmul-error from empty matrix 
	if (m == 0) {
		std::cout << "INFO: No voxel found in [ " << fn << " ], program terminating .." << std::endl;
		exit(0);
	}
	double* q = new double[nq];
	double* intensity = new double[nq];
	for (i = 0; i < nq; i++) q[i] = 0.002 * i;

	// Solid angle average
	double* xx = new double[J];
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

	//array mat(nr_voxels, 3, f64);
	//mat = coordata; 
	array qr(m, J, f32);
	array sys1(1, J, f32), sys2(1, J, f32);

	std::cout << "Start to calculate swaxs curves..." << std::endl;
	timer::start();

	// Matrix operations using GPU for every q[i]
	for (i = 0; i < nq; i++) {
		qr = q[i] * matmul(coordata, qmat);
		sys1 = sum(cos(qr));
		// NOTE: Here, a minus sign does not matter, squared out anyway
		sys2 = sum(1.0 * sin(qr));
		sys2 = pow(sys2, 2) + pow(sys1, 2);
		intensity[i] = mean(sys2, 1)(0).scalar<float>();

		// Using A Cool Loading Bar
		std::cout << "\r" << (100 * i / (nq - 1)) << "% completed: ";
		std::cout << std::string((50 * i / nq), '|');
		std::cout.flush();
	}

	printf("\nINFO: Elapsed time: %g seconds\n", timer::stop());

	// Print out the data file
	// The intensities are normalized using I(0) with 17 digits precision
	filename.replace(filename.end() - 6, filename.end(), "dat");
	std::ofstream* out = new std::ofstream(filename);
	if (!out->good()) {
		std::cout << "Error: Cannot open [ " << filename << " ]. " << std::endl;
		exit(1);
	}

	std::cout << "INFO: Writing q, and intensity to the file: " << filename << " ." << std::endl;
	(*out).precision(17);
	for (int i = 0; i < nq; i++) {
		*out << q[i] << "\t" << intensity[i] / intensity[0] << std::endl;
	}
	out->close();
	std::cout << "INFO: Writing completed successfully." << std::endl;

	delete[] q;
	delete[] intensity;
}