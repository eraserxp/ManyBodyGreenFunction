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
	nmax;
	t0;
	d0;
	std::vector<int> impurityIndex; // the index for the impurity site
	std::vector<double> impurityStrength; //disorder_strength;
};

void formBasisSets(int nmax, IMatrix& indexMatrix, PairVector& basisSets);
void formHamiltonianMatrix(InputParameters& pars, DMatrix& hamiltonian, IMatrix& indexMatrix,
		                   PairVector& basisSets);
void obtainEigenVectors(DMatrix& hamiltonian, DVector& eigenValues, DMatrix& eigenVectors);
#endif /* DIRECT_CALCULATION_H_ */
