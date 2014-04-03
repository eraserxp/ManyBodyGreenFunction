/*
 * scattering.cpp
 *
 *  Created on: Jan 29, 2014
 *      Author: pxiang
 */

#include "scattering.h"

void printMatrix(CDMatrix& m) {
	int rows = m.rows();
	int cols = m.cols();
	std::cout << "The content of the matrix is:" << std::endl;
	for (int r=0; r<rows; ++r) {
		for (int c=0; c<cols; ++c) {
			std::cout << r << "  " << c << "  "<< m(r,c) <<std::endl;
		}
	}
}

void printVector(CDVector& v) {
	int rows = v.rows();
	std::cout << "The content of the vector is:" << std::endl;
	for (int r=0; r<rows; ++r) {
		std::cout << r << "  " << v(r) <<std::endl;

	}
}

void normalizeWavefunction(CDVector& v) {
	DVector square = v.cwiseAbs2();
	double sum = square.sum();
	double normalizationFactor = 1.0/sqrt(sum);
	v = v*dcomplex(normalizationFactor, 0.0);
}

void changeElements(CDMatrix& cm) {
	int rows = cm.rows();
	int cols = rows;
//	#pragma omp parallel for
	for (int r=0; r<rows; ++r) {
		for (int c=0; c<cols; ++c) {
			cm(r,c) = -cm(r,c);
		}
	}

//	#pragma omp parallel for
	for (int i=0; i<rows; ++i) {
		cm(i,i) += 1.0;
	}
}

double biexcitonEnergy(double K, double t, double d) {
	d = 2.0*d;
	return d + 4*t*t*cos(K/2.0)*cos(K/2.0)/d;
}

// the wavefunction of biexciton is given by \sum_{n, m}exp[i*K*(n+m)/2]*psi(|m-n|) |n, m>
double psi(int separation, double K, double t, double d, int nmol) {
	double cosine = cos(K/2.0);
	d = 2.0*d;
	return sqrt(d*d - 4*t*t*cosine*cosine)/(2.0*d*sqrt(nmol*1.0))*
			pow(2*t*cosine/d, separation-1);
}


// calculate the initial state which is a bound pair of excitations
// record the basis set (i,j) in biexcitonBasis
void generateInitialState(Parameters& pars,  int impuritySite, double K, PairVector& biexcitonBasis,
		    PairVector& disorderBasis, ComplexMatrix& initialState,
		    IntegerMatrix& disorderBasisPosition) {

	dcomplex ii = dcomplex(0.0, 1.0);
	int nmax = pars.nmax;

	// calculate separation_max
	int separation_max = 1;
	double base = 2*pars.t0*cos(K/2.0)/pars.d0;
	double decayFactor = base;
	do {
		decayFactor = decayFactor*base;
		separation_max++;
	} while (decayFactor<base*0.01);

	int rows = (nmax - separation_max + 1)*separation_max; // the dimension of biexciton state
	initialState = ComplexMatrix(rows, 1);
	// to record the position of disorder basis in biexciton basis
	disorderBasisPosition = IntegerMatrix(nmax, separation_max+1);

	int nth = 0;
	biexcitonBasis.clear();
	disorderBasis.clear();

	for (int i=0; i<=nmax-separation_max; ++i) {
		for (int sep=1; sep<=separation_max; ++sep) {
			int j = i + sep;
			Pair basis(i,j);
			biexcitonBasis.push_back(basis);
			if (i==impuritySite || j==impuritySite) { // if (i,j) belongs to disorderBasis
				disorderBasis.push_back(basis);
				disorderBasisPosition(i,sep) = nth;
			}
			dcomplex tmp = std::exp(ii*(i*1.0+j*1.0)/2.0)* \
					psi(sep, K, pars.t0, pars.d0, pars.nmax + 1);
			initialState(nth,0).real = tmp.real();
			initialState(nth,0).imag = tmp.imag();
			nth++;
		}
	}

}



// calculate the initial state which is a bound pair of excitations
// record the basis set (i,j) in biexcitonBasis
void generateInitialState2(Parameters& pars,  int impuritySite, double K, PairVector& biexcitonBasis,
		    PairVector& disorderBasis, CDVector& initialState,
		    IMatrix& disorderBasisPosition) {

	dcomplex ii = dcomplex(0.0, 1.0);
	int nmax = pars.nmax;

	// calculate separation_max
	int separation_max = 1;
	double d0 = pars.d0*2;
	double base = 2*pars.t0*cos(K/2.0)/d0;
	double decayFactor = base;
	do {
		decayFactor = decayFactor*base;
		separation_max++;
	} while (decayFactor<base*0.01);

	int rows = (nmax - separation_max + 1)*separation_max; // the dimension of biexciton state
	initialState = CDVector(rows);
	// to record the position of disorder basis in biexciton basis
	disorderBasisPosition.resize(nmax, separation_max+1);

	int nth = 0;
	biexcitonBasis.clear();
	disorderBasis.clear();

	for (int i=0; i<=nmax-separation_max; ++i) {
		for (int sep=1; sep<=separation_max; ++sep) {
			int j = i + sep;
			Pair basis(i,j);
			biexcitonBasis.push_back(basis);
			if (i==impuritySite || j==impuritySite) { // if (i,j) belongs to disorderBasis
				disorderBasis.push_back(basis);
				disorderBasisPosition(i,sep) = nth;
			}
			dcomplex tmp = 2.0*std::exp(ii*K*(i*1.0+j*1.0)/2.0)* \
					psi(sep, K, pars.t0, pars.d0, pars.nmax + 1);
			initialState(nth) = tmp;
			nth++;
		}
	}

}

//
void calculateScatteringState(Parameters& pars, double K, int impuritySite, double disorderStrength,
		                      ComplexMatrix& scatteringState, double& transmissionCoeff) {

	PairVector biexcitonBasis;
	PairVector disorderBasis;
	ComplexMatrix initialState;
	IntegerMatrix disorderBasisPosition;
//	printf("1.1\n");

	// calculate the initial state
	generateInitialState(pars,  impuritySite, K,  biexcitonBasis,
			    disorderBasis, initialState, disorderBasisPosition);

//	printf("1.2\n");
	AlphaBeta ab(pars);
	IntegerMatrix indexMatrix;
	ab.FillIndexMatrix(indexMatrix);

//	printf("1.3\n");

	int biexciton_size = biexcitonBasis.size();
	int disorder_size = disorderBasis.size();
	//green's function in the biexciton basis and disorder basis
	ComplexMatrix gf_biexciton(biexciton_size, disorder_size);

	//green's function in the disorder basis
	ComplexMatrix gf_disorder(disorder_size, disorder_size);

	ComplexMatrix h1(disorder_size, disorder_size);
	complex_mkl z;
	z.real = biexcitonEnergy(K,pars.t0,pars.d0); // the biexciton energy is determined by K
	z.imag = 0.01;

//	printf("1.4\n");
	for (int c=0; c<disorder_size; ++c) {
		Pair basis = disorderBasis[c];
		int ni1 = basis.First();
		int ni2 = basis.Second();
//		printf("1.5\n");
		// calculate all elements of Green's function and save them into disk
		// ni1 and ni2 are the initial positions of the two excitations
		calculateAllGF(ni1, ni2,z,ab);
//		printf("1.6\n");

		for (int r=0; r<biexciton_size; ++r) {
			basis = biexcitonBasis[r];
			int nf1 = basis.First();
			int nf2 = basis.Second();
			// calculate the green's function in the biexciton basis
			gf_biexciton(r,c) = extractMatrixElement(nf1, nf2, ni1, ni2, indexMatrix);
		}

//		printf("1.7\n");

		for (int r=0; r<disorder_size; ++r) {
			// calculate the green's function in the disorder basis
			basis = disorderBasis[r];
			int nf1 = basis.First();
			int nf2 = basis.Second();
			gf_disorder(r,c) = extractMatrixElement(nf1, nf2, ni1, ni2, indexMatrix);
			// form matrix h1
			if (r==c) {
				h1(r,c).real = disorderStrength;
				h1(r,c).imag = 0.0;
			} else {
				h1(r,c).real = 0.0;
				h1(r,c).imag = 0.0;
			}
		}
//		printf("1.8\n");
	}

//	printf("1.9\n");
	// T matrix in disorder basis
	ComplexMatrix TMatrix = gf_disorder*h1;
//	printf("1.10\n");
	changeElements(TMatrix);//calculate 1 - gf_disorder*h1
//	printf("1.11\n");
	TMatrix.InverseInPlace();
//	printf("1.12\n");
	ComplexMatrix tmp2 = h1*TMatrix;
//	printf("1.13\n");
	TMatrix = tmp2;

	ComplexMatrix initial_disorder(disorder_size, 1); // initial state in disorder basis
//	printf("1.14\n");
	for (int i=0; i<disorder_size; ++i) {
		Pair basis = disorderBasis[i];
		int n = basis.First();
		int m = basis.Second();
		int nth = disorderBasisPosition(n, m-n);
		initial_disorder(i,0) = initialState(nth, 0);
	}
//	printf("1.15\n");
	ComplexMatrix beta_T = TMatrix*initial_disorder;
//	printf("1.16\n");
	ComplexMatrix beta_GT = gf_biexciton*beta_T;
//	printf("1.17\n");
	ComplexMatrix tmp = initialState + beta_GT;
	//scatteringState = initialState + beta_GT;
//	printf("1.18\n");
	scatteringState = tmp;
//	printf("1.19\n");

	// calculate the transmission coefficient at the position far away from the disorder site
	int nth = disorderBasisPosition(pars.nmax-6,1);
	dcomplex  original = dcomplex(initialState(nth,0).real, initialState(nth,0).imag);
	dcomplex  transmitted = dcomplex(scatteringState(nth,0).real, scatteringState(nth,0).imag);
	dcomplex ratio = transmitted/original;
	transmissionCoeff = std::abs(ratio)*std::abs(ratio);
}



void calculateScatteringState2(Parameters& pars, double K, int impuritySite, double disorderStrength,
		                      CDVector& scatteringState, double& transmissionCoeff) {

	PairVector biexcitonBasis;
	PairVector disorderBasis;
	CDVector initialState;
	IMatrix disorderBasisPosition;
//	printf("1.1\n");

	// calculate the initial state
	generateInitialState2(pars,  impuritySite, K,  biexcitonBasis,
			    disorderBasis, initialState, disorderBasisPosition);

//	printf("1.2\n");
	AlphaBeta ab(pars);
	IntegerMatrix indexMatrix;
	ab.FillIndexMatrix(indexMatrix);

//	printf("1.3\n");

	int biexciton_size = biexcitonBasis.size();
	int disorder_size = disorderBasis.size();
	//green's function in the biexciton basis and disorder basis
	CDMatrix gf_biexciton(biexciton_size, disorder_size);

	//green's function in the disorder basis
	CDMatrix gf_disorder(disorder_size, disorder_size);

	CDMatrix h1(disorder_size, disorder_size);
	complex_mkl z;
	z.real = biexcitonEnergy(K,pars.t0,pars.d0); // the biexciton energy is determined by K
	z.imag = 0.0001;

//	printf("1.4\n");
	for (int c=0; c<disorder_size; ++c) {
		Pair basis = disorderBasis[c];
		int ni1 = basis.First();
		int ni2 = basis.Second();
//		printf("1.5\n");
		// calculate all elements of Green's function and save them into disk
		// ni1 and ni2 are the initial positions of the two excitations
		calculateAllGF(ni1, ni2,z,ab);
//		printf("1.6\n");

		for (int r=0; r<biexciton_size; ++r) {
			basis = biexcitonBasis[r];
			int nf1 = basis.First();
			int nf2 = basis.Second();
			// calculate the green's function in the biexciton basis
			complex_mkl tmp = extractMatrixElement(nf1, nf2, ni1, ni2, indexMatrix);
			gf_biexciton(r,c) = dcomplex(tmp.real, tmp.imag);

		}

//		printf("1.7\n");

		for (int r=0; r<disorder_size; ++r) {
			// calculate the green's function in the disorder basis
			basis = disorderBasis[r];
			int nf1 = basis.First();
			int nf2 = basis.Second();
			complex_mkl tmp = extractMatrixElement(nf1, nf2, ni1, ni2, indexMatrix);
			gf_disorder(r,c) = dcomplex(tmp.real, tmp.imag);
			//std::cout << r << "  " << c << "  " << gf_disorder(r,c) << std::endl;
			// form matrix h1
			if (r==c) {
				h1(r,c) = dcomplex(disorderStrength, 0.0);
				//std::cout << "h1 " << r << "  " << c << "  " << h1(r,c) << std::endl;
			} else {
				h1(r,c) = dcomplex(0.0, 0.0);

			}
		}
//		printf("1.8\n");
	}

//	printf("1.9\n");
	// T matrix in disorder basis
	CDMatrix TMatrix = gf_disorder*h1;
//	printf("1.10\n");
	changeElements(TMatrix);//calculate 1 - gf_disorder*h1
//	printf("1.11\n");
	TMatrix = TMatrix.inverse();
//	printf("1.12\n");
	CDMatrix tmp2 = h1*TMatrix;
//	printf("1.13\n");
	TMatrix = tmp2;

	//printMatrix(TMatrix);

	CDVector initial_disorder(disorder_size); // initial state in disorder basis
	// normalize initial wavefunction
	normalizeWavefunction(initialState);

//	printf("1.14\n");
	for (int i=0; i<disorder_size; ++i) {
		Pair basis = disorderBasis[i];
		int n = basis.First();
		int m = basis.Second();
		int nth = disorderBasisPosition(n, m-n);
		initial_disorder(i) = initialState(nth);
	}
//	printf("1.15\n");
	CDVector beta_T = TMatrix*initial_disorder;

	//printVector(beta_T);

//	printf("1.16\n");
	CDVector beta_GT = gf_biexciton*beta_T;
//	printf("1.17\n");
	CDVector tmp = initialState + beta_GT;
	//scatteringState = initialState + beta_GT;
//	printf("1.18\n");
	scatteringState = tmp;
//	printf("1.19\n");

	// normalize wavefunction
	normalizeWavefunction(scatteringState);

	// calculate the transmission coefficient at the position far away from the disorder site
	int nth = disorderBasisPosition(pars.nmax-8,1);
	dcomplex  original = initialState(nth);
	dcomplex  transmitted = scatteringState(nth);
	dcomplex ratio = transmitted/original;
	transmissionCoeff = std::abs(ratio)*std::abs(ratio);
}



// calculate the scattering state making use of direct diagonalization
//void calculateScatteringState_direct(InputParameters& pars, DVector& transmissionCoeff) {
//
//	int nmax = pars.nmax;
//	IMatrix indexMatrix;
//	PairVector basisSets;
//	formBasisSets(nmax, indexMatrix, basisSets);
//
//	DMatrix hamiltonian;
//	formHamiltonianMatrix(pars, hamiltonian, indexMatrix, basisSets);
//
//	DVector eigenValues;
//	DMatrix eigenVectors;
//	obtainEigenVectors(hamiltonian, eigenValues, eigenVectors);
//
//	int n1f = n1i;
//	int n2f = n2i;
//	CDArray numerator;
//	numeratorHelper(n1f, n2f, n1i, n2i, indexMatrix, eigenVectors, numerator);
//
//	dosList.clear();
//	dosList.reserve(zList.size());
//	for (int i=0; i<zList.size(); ++i) {
//		dcomplex z = zList[i];
//		CDArray oneOverDenominator;
//		denominatorHelper(z, eigenValues, oneOverDenominator);
//		dcomplex gf = greenFunc(numerator, oneOverDenominator);
//		dosList[i] = -gf.imag()/M_PI;
//	}
//
//	// calculate the transmission coefficient at the position far away from the disorder site
//	int nth = disorderBasisPosition(pars.nmax-8,1);
//	dcomplex  original = initialState(nth);
//	dcomplex  transmitted = scatteringState(nth);
//	dcomplex ratio = transmitted/original;
//	transmissionCoeff = std::abs(ratio)*std::abs(ratio);
//}

