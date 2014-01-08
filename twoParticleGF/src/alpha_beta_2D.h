/*
 * generate_alpha_beta_2D.h
 *
 *  Created on: Jan 6, 2014
 *      Author: pxiang
 */
#ifndef ALPHA_BETA_2D_H_
#define ALPHA_BETA_2D_H_

#include "types.h"
#include "mics.h"
#include "complex_matrix.h"
#include "integer_matrix.h"

#include "random_generator.h"
#include "alpha_beta.h"


#include <vector>


// generate the hamiltonian matrix
class Hamiltonian2D {
public:

	Hamiltonian2D(Parameters2D& p);

	double energyAtSite(int nx, int ny);

	double t(int nx, int ny, char direction);

	double d(int nx, int ny, char direction);

private:
	int xmax;
	int ymax;
	double e0, t0, d0;
	double e0MaxDisorder, t0MaxDisorder, d0MaxDisorder;
	unsigned e0seed, t0seed, d0seed;
	OnsiteDisorder  rnm_e0;
	RandomInteraction2D  rnm_t0;
	RandomInteraction2D  rnm_d0;

};


void getCoordinates(int index, int xmax, int ymax, int& xi, int& yi);

int coordinatesToIndex(int xmax, int ymax, int xi, int yi);

Neighbors2D generateNeighbors2D(int xi, int yi, int xj, int yj, int xmax, int ymax);

void generateIndexMatrix2D(int xmax, int ymax, IntegerMatrix& Index, QuartetListVector& VtoG,
		std::vector<int>& DimsOfV);

class AlphaBeta2D {
public:
	AlphaBeta2D(Parameters2D& ps);
	void FillAlphaBetaMatrix(int nsum, complex_mkl z, ComplexMatrix& alpha, ComplexMatrix& beta);
	int GetXmax();
	int GetYmax();
	int GetIndex(int x1, int y1, int x2, int y2);
	int GetDimOfV(int nsum);
	// obtain the factor in front of G(ni1, ni2) in the equations for Green's functions
	complex_mkl GetFactor(int x1_i, int y1_i, int x2_i, int y2_i, complex_mkl z);

private:
	int xmax;
	int ymax;
	double e0, t0, d0;
	double e0MaxDisorder, t0MaxDisorder, d0MaxDisorder;
	unsigned e0seed, t0seed, d0seed;
	Hamiltonian2D ham;
	IntegerMatrix Index;
	QuartetListVector VtoG;
	std::vector<int> DimsOfV;
};

ComplexMatrix fromRightToCenter2D(int Kc, complex_mkl z, AlphaBeta2D& ab);
ComplexMatrix fromLeftToCenter2D(int Kc, complex_mkl z, AlphaBeta2D& ab);
ComplexMatrix solveVnc2D(int x1_i, int y1_i, int x2_i, int y2_i, complex_mkl z, AlphaBeta2D& ab);

void generateDensityOfState2D(int x1_i, int y1_i, int x2_i, int y2_i, Parameters2D& pars,
		        const std::vector<complex_mkl>& zList, std::vector<double>& rhoList);

#endif
