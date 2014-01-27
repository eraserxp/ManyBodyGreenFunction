/*
 * alpha_beta_sparse.h
 *
 *  Created on: Jan 27, 2014
 *      Author: pxiang
 */

#ifndef ALPHA_BETA_SPARSE_H_
#define ALPHA_BETA_SPARSE_H_

#include "types.h"
#include "mics.h"
#include "complex_matrix.h"
#include "integer_matrix.h"

#include "random_generator.h"
#include "alpha_beta.h"


#include <vector>

// use sparse matrix from eigen c++ library to do the calculations
class AlphaBetaSP {
public:
	AlphaBetaSP(Parameters& ps);
	void FillAlphaBetaMatrix(int nsum, dcomplex z, CSMatrix& alpha, CSMatrix& beta);
	int GetNmax();
	// obtain the factor in front of G(ni1, ni2) in the equations for Green's functions
	dcomplex GetFactor(int ni1, int ni2, dcomplex z);

private:
	int nmax;
	double e0, t0, d0;
	double e0MaxDisorder, t0MaxDisorder, d0MaxDisorder;
	unsigned e0seed, t0seed, d0seed;
	Hamiltonian ham;
	IntegerMatrix Index;
	PairMatrix VtoG;
	std::vector<int> DimsOfV;
};


void changeElements(CSMatrix& csm);

CDMatrix solveDenseLinearEqs(CDMatrix& A, CDMatrix& B);

CSMatrix solveDenseLinearEqs(CSMatrix& A, CSMatrix& B);

CSMatrix solveSparseLinearEqs(CSMatrix& A, CSMatrix& B);

CDVector solveSparseLinearEqs2(CSMatrix& A, CDVector& b);


CSMatrix fromRightToCenterSP(int Kc, dcomplex z, AlphaBetaSP& ab);
void checkSparsenessRightToCenter(int Kc, dcomplex z, AlphaBetaSP& ab);

CSMatrix fromLeftToCenterSP(int Kc, dcomplex z, AlphaBetaSP& ab);
void checkSparsenessLeftToCenter(int Kc, dcomplex z, AlphaBetaSP& ab);

void checkSparseness();


//
//
CDVector solveVncSP(int ni1, int ni2, dcomplex z, AlphaBetaSP& ab);
//
void generateDensityOfStateSP(int ni1, int ni2, Parameters& pars, const std::vector<dcomplex>& zList, std::vector<double>& rhoList);


#endif /* ALPHA_BETA_SPARSE_H_ */
