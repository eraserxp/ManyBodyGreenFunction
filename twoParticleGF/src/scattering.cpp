/*
 * scattering.cpp
 *
 *  Created on: Jan 29, 2014
 *      Author: pxiang
 */

#include "scattering.h"

double biexcitonEnergy(double K, Parameters pars) {
	return pars.d0 + 4*pars.t0*pars.t0*cos(K/2.0)*cos(K/2.0)/pars.d0;
}

// the wavefunction of biexciton is given by \sum_{n, m}exp[i*K*(n+m)/2]*psi(|m-n|) |n, m>
double psi(int separation, double K, double t, double d, int nmol) {
	double cosine = cos(K/2.0);
	return sqrt(d*d - 4*t*t*cosine*cosine)/(2.0*d*sqrt(nmol*1.0))*
			pow(2*t*cosine/d, separation-1);
}


// calculate the initial state which is a bound pair of excitations
// record the basis set (i,j) in biexcitonBasis
void generateInitialState(int nmax,  int impuritySite, double K,
		Parameters& pars, PairVector& biexcitonBasis, ComplexMatrix& initialState) {
	dcomplex ii = dcomplex(0.0, 1.0);
	int separation_max = 1;
	// calculate separation_max
	double base = 2*pars.t0*cos(K/2.0)/pars.d0;
	double decayFactor = base;
	do {
		decayFactor = decayFactor*base;
		separation_max++;
	} while (decayFactor<base*0.01);

	int rows = (nmax - separation_max + 1)*separation_max;
	initialState = ComplexMatrix(rows, 1);
	int index = 0;
	biexcitonBasis.clear();
	for (int i=0; i<=nmax-separation_max; ++i) {
		for (int sep=1; sep<=separation_max; ++sep) {
			int j = i + sep;
			Pair basis(i,j);
			biexcitonBasis.push_back(basis);
			dcomplex tmp = std::exp(ii*(i*1.0+j*1.0)/2.0)* \
					psi(sep, K, pars.t0, pars.d0, pars.nmax + 1);
			initialState(index,0).real = tmp.real();
			initialState(index,0).imag = tmp.imag();
			index++;
		}
	}

}


// exchange row1 and row2 of a complex matrix
void exchangeTwoRows(ComplexMatrix& m, int row1, int row2) {
	int cols = m.GetCols();
	complex_mkl tmp;

	for (int c=0; c<cols; ++c) {
		tmp = m(row1, c);
		m(row1, c) = m(row2, c);
		m(row2, c) = tmp;
	}
}
