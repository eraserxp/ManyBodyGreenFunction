/*
 * direct_calculation_test.cpp
 *
 *  Created on: Jan 12, 2014
 *      Author: pxiang
 */

#include "direct_calculation.h"
#include "gtest/gtest.h"
#include "mics.h"

TEST(DirectCalculationTest, DISABLED_CheckDOS) {
	InputParameters pars;
	pars.nmax = 100;
	pars.t0 = 5.0;
	pars.d0 = 15.0;
	pars.impurityIndex.clear();
	pars.impurityStrength.clear();

	int n1i = pars.nmax/2;
	int n2i = n1i + 1;
	std::vector<dcomplex > zList(101);
	std::vector<double> zRealList = linspace(-50,50,101);

	for (int i=0; i<zList.size(); ++i) {
		zList[i] =dcomplex( zRealList[i], 0.1);
	}

	std::vector<double> rhoList;
	densityOfState_direct(pars, n1i, n2i, zList, rhoList);
	save_two_arrays("rho_vs_energy_direct.txt", zRealList, rhoList);
	EXPECT_EQ(2,2);
}


TEST(DirectCalculationTest, DISABLED_CheckDOSTwoHundred) {
	InputParameters pars;
	pars.nmax = 200;
	pars.t0 = 5.0;
	pars.d0 = 15.0;
	pars.impurityIndex.clear();
	pars.impurityStrength.clear();

	int n1i = pars.nmax/2;
	int n2i = n1i + 1;
	std::vector<dcomplex > zList(101);
	std::vector<double> zRealList = linspace(-50,50,101);

	for (int i=0; i<zList.size(); ++i) {
		zList[i] =dcomplex( zRealList[i], 0.1);
	}

	std::vector<double> rhoList;
	densityOfState_direct(pars, n1i, n2i, zList, rhoList);
	save_two_arrays("rho_vs_energy_direct_200.txt", zRealList, rhoList);
	EXPECT_EQ(2,2);
}
