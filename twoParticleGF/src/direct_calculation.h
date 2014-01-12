/*
 * direct_calculation.h
 *
 *  Created on: Jan 10, 2014
 *      Author: pxiang
 */

#ifndef DIRECT_CALCULATION_H_
#define DIRECT_CALCULATION_H_

#include "types.h"
#include "alpha_beta.h"

struct InputParameters {
	int nmax;
	double t0;
	double d0;
	std::vector<int> impurityIndex; // the index for the impurity site
	std::vector<double> impurityStrength; //disorder_strength;
};

void formBasisSets(int nmax, IMatrix& indexMatrix, PairVector& basisSets);
void formHamiltonianMatrix(InputParameters& pars, DMatrix& hamiltonian, IMatrix& indexMatrix,
		                   PairVector& basisSets);
void obtainEigenVectors(DMatrix& hamiltonian, DVector& eigenValues, DMatrix& eigenVectors);

void numeratorHelper(int n1f, int n2f, int n1i, int n2i, IMatrix& indexMatrix,
		             DMatrix& eigenVectors, CDArray& numerator);

void denominatorHelper(dcomplex z, DVector& eigenValues, CDArray& denominator);

dcomplex greenFunc(CDArray& numerator, CDArray& oneOverDenominator);

void densityOfState_direct(InputParameters& pars, int ni1, int ni2, std::vector<dcomplex >& zList,
		                   std::vector<double>& dosList);
#endif /* DIRECT_CALCULATION_H_ */
