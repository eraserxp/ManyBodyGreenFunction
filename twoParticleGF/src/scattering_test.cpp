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
	pars.t0 = 5.0;
	pars.d0 = 15.0;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;

	std::vector<double> K_array = linspace(-M_PI,M_PI,2);
	std::vector<double> transmission_array;
	int impuritySite = pars.nmax/2;
	double disorderStrength = 5.0;
	ComplexMatrix scatteringState;
	double transmission;
	for (int i=0; i<K_array.size();++i) {
		double K = K_array[i];
		calculateScatteringState(pars, K, impuritySite, disorderStrength, scatteringState, transmission);
		transmission_array.push_back(transmission);
	}
	// save to the disk
	save_two_arrays("transmission_vs_k.txt",K_array, transmission_array);

	EXPECT_EQ(2,2);
}
