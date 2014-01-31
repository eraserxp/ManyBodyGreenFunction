/*
 * scattering.h
 *
 *  Created on: Jan 29, 2014
 *      Author: pxiang
 */

#ifndef SCATTERING_H_
#define SCATTERING_H_

#include "types.h"
#include <math.h>       /* cos */
#include "complex_matrix.h"
#include "integer_matrix.h"
#include "alpha_beta.h"

double biexcitonEnergy(double K, double t, double d);

double psi(int separation, double K, double t, double d, int nmol);

void generateInitialState(Parameters& pars,  int impuritySite, double K, PairVector& biexcitonBasis,
		    PairVector& disorderBasis, ComplexMatrix& initialState,
		    IntegerMatrix& disorderBasisPosition);

void calculateScatteringState(Parameters& pars, double K, int impuritySite, double disorderStrength,
		                      ComplexMatrix& scatteringState, double& transmissionCoeff);

#endif /* SCATTERING_H_ */
