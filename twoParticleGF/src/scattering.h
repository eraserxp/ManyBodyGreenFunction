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

void printMatrix(CDMatrix& m);

void normalizeWavefunction(CDVector& v);

void changeElements(CDMatrix& cm);

double biexcitonEnergy(double K, double t, double d);

double psi(int separation, double K, double t, double d, int nmol);

// use my own matrix class (there are memory leaks somewhere)
void generateInitialState(Parameters& pars,  int impuritySite, double K, PairVector& biexcitonBasis,
		    PairVector& disorderBasis, ComplexMatrix& initialState,
		    IntegerMatrix& disorderBasisPosition);

void calculateScatteringState(Parameters& pars, double K, int impuritySite, double disorderStrength,
		                      ComplexMatrix& scatteringState, double& transmissionCoeff);


// use the matrix class from eigen c++
void generateInitialState2(Parameters& pars,  int impuritySite, double K, PairVector& biexcitonBasis,
		    PairVector& disorderBasis, CDVector& initialState,
		    IMatrix& disorderBasisPosition);

void calculateScatteringState2(Parameters& pars, double K, int impuritySite, double disorderStrength,
		                      CDVector& scatteringState, double& transmissionCoeff);

#endif /* SCATTERING_H_ */
