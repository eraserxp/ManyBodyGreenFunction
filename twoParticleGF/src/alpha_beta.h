/*
 * generate_alpha_beta.h
 *
 *  Created on: Dec 18, 2013
 *      Author: pxiang
 */

#ifndef GENERATE_ALPHA_BETA_H_
#define GENERATE_ALPHA_BETA_H_

#include "types.h"
#include "mics.h"
#include "complex_matrix.h"
#include "integer_matrix.h"

#include "random_generator.h"


#include <vector>

dcomplex convertToDcomplex(complex_mkl& mkl_complex);

// generate the hamiltonian matrix
class Hamiltonian {
public:

	Hamiltonian(Parameters& p);

	double energyAtSite(int n);

	double t(int n, int m);

	double d(int n, int m);

private:
	int nmax;
	double e0, t0, d0;
	double e0MaxDisorder, t0MaxDisorder, d0MaxDisorder;
	unsigned e0seed, t0seed, d0seed;
	RandomNumberMatrix  rnm_e0;
	RandomNumberMatrix  rnm_t0;
	RandomNumberMatrix  rnm_d0;

};





void generateIndexMatrix(int nmax, IntegerMatrix& Index, PairMatrix& VtoG, std::vector<int>& DimsOfV);

Neighbors generateNeighbors(int i, int j, int nmax);

class AlphaBeta {
public:
	AlphaBeta(Parameters& ps);
	void FillAlphaBetaMatrix(int nsum, complex_mkl z, ComplexMatrix& alpha, ComplexMatrix& beta);
	int GetNmax();
	// obtain the factor in front of G(ni1, ni2) in the equations for Green's functions
	complex_mkl GetFactor(int ni1, int ni2, complex_mkl z);

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

void changeElements(ComplexMatrix& cm);

ComplexMatrix fromRightToCenter(int Kc, complex_mkl z, AlphaBeta& ab);
ComplexMatrix fromLeftToCenter(int Kc, complex_mkl z, AlphaBeta& ab);


ComplexMatrix solveVnc(int ni1, int ni2, complex_mkl z, AlphaBeta& ab);

// given g(ni,nj), it tells the index of g(ni,nj) in vector V_{ni+nj}
// if g(ni,nj) is not in V_{ni+nj}, it returns -1 which means not found
int getIndex(int nmax, int nsum, int ni, int nj);


void generateDensityOfState(int ni1, int ni2, Parameters& pars, const std::vector<complex_mkl>& zList, std::vector<double>& rhoList);






















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

double dos(double E, void* pars);

double func_test(double E, void* pars);

#endif /* GENERATE_ALPHA_BETA_H_ */
