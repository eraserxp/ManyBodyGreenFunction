/*
 * scattering_test.cpp
 *
 *  Created on: Jan 31, 2014
 *      Author: pxiang
 */

#include "scattering.h"
#include "gtest/gtest.h"
#include "mics.h"
#include "types.h"

TEST(ScatteringStateTest, OneHundred) {
	Parameters pars;
	pars.nmax = 100;
	pars.e0 = 0.0;
	pars.t0 = 1.0;
	pars.d0 = 3.0;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;

	std::vector<double> K_array = linspace(-M_PI*0.99,0.0,10);
	std::vector<double> transmission_array;
	int impuritySite = pars.nmax/2;
	double disorderStrength = 5.0;


	for (int i=0; i<K_array.size();++i) {
		ComplexMatrix scatteringState;
		double transmission;
		double K = K_array[i];
		calculateScatteringState(pars, K, impuritySite, disorderStrength, scatteringState, transmission);
		transmission_array.push_back(transmission);
	}
	// save to the disk
	save_two_arrays("transmission_vs_k.txt",K_array, transmission_array);

	EXPECT_EQ(2,2);
}


TEST(GenerateInitialState, OneHundred) {
	Parameters pars;
	pars.nmax = 100;
	pars.e0 = 0.0;
	pars.t0 = 1.0;
	pars.d0 = 3.0;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;
	int impuritySite = pars.nmax/2;

	PairVector biexcitonBasis;
	PairVector disorderBasis;
	ComplexMatrix initialState;
	IntegerMatrix disorderBasisPosition;
	//printf("1.1\n");

	std::vector<double> K_array = linspace(-M_PI*0.99,M_PI*0.99,100);
	for (int i=0; i<K_array.size();++i) {
		double K = K_array[i];
		// calculate the initial state
		generateInitialState(pars,  impuritySite, K,  biexcitonBasis,
			    disorderBasis, initialState, disorderBasisPosition);

	}

	EXPECT_EQ(2,2);
}
