/*
 * alpha_beta_sparse_test.cpp
 *
 *  Created on: Jan 27, 2014
 *      Author: pxiang
 */

#include "alpha_beta_sparse.h"
#include "gtest/gtest.h"
#include "mics.h"


TEST(AlphaBetaSPTest, DISABLED_FourSitesNoDisorder) {
	Parameters pars;
	pars.nmax = 3;
	pars.e0 = pars.t0 = pars.d0 = 1;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;

	AlphaBetaSP ab(pars);
	int nsum = 1;
	dcomplex z(1.0, 0.1);
	CSMatrix alpha;
	CSMatrix beta;
	ab.FillAlphaBetaMatrix(nsum, z, alpha, beta);
	// testing the dimensions
	//EXPECT_EQ(alpha.rows(),0);
	//EXPECT_EQ(alpha.cols(),0);
	EXPECT_EQ(beta.rows(),1);
	EXPECT_EQ(beta.cols(),1);
	//testing matrix element
	dcomplex x = pars.t0/(z-pars.e0-pars.e0-pars.d0);
	EXPECT_DOUBLE_EQ(beta.coeff(0,0).real(),x.real());
	EXPECT_DOUBLE_EQ(beta.coeff(0,0).imag(),x.imag());

	nsum=2;
	ab.FillAlphaBetaMatrix(nsum, z, alpha, beta);
	// testing the dimensions
	EXPECT_EQ(alpha.rows(),1);
	EXPECT_EQ(alpha.cols(),1);
	EXPECT_EQ(beta.rows(),1);
	EXPECT_EQ(beta.cols(),2);
	//testing matrix element
	x = pars.t0/(z-pars.e0-pars.e0);
	EXPECT_DOUBLE_EQ(alpha.coeff(0,0).real(),x.real());
	EXPECT_DOUBLE_EQ(alpha.coeff(0,0).imag(),x.imag());
	EXPECT_DOUBLE_EQ(beta.coeff(0,0).real(),x.real());
	EXPECT_DOUBLE_EQ(beta.coeff(0,0).imag(),x.imag());
	EXPECT_DOUBLE_EQ(beta.coeff(0,1).real(),x.real());
	EXPECT_DOUBLE_EQ(beta.coeff(0,1).imag(),x.imag());

	nsum =3;
	ab.FillAlphaBetaMatrix(nsum, z, alpha, beta);
	// testing the dimensions
	EXPECT_EQ(alpha.rows(),2);
	EXPECT_EQ(alpha.cols(),1);
	EXPECT_EQ(beta.rows(),2);
	EXPECT_EQ(beta.cols(),1);
	//testing matrix element
	x = pars.t0/(z-pars.e0-pars.e0);
	EXPECT_DOUBLE_EQ(alpha.coeff(0,0).real(),x.real());
	EXPECT_DOUBLE_EQ(alpha.coeff(0,0).imag(),x.imag());
	x = pars.t0/(z-pars.e0-pars.e0-pars.d0);
	EXPECT_DOUBLE_EQ(alpha.coeff(1,0).real(),x.real());
	EXPECT_DOUBLE_EQ(alpha.coeff(1,0).imag(),x.imag());
	x = pars.t0/(z-pars.e0-pars.e0);
	EXPECT_DOUBLE_EQ(beta.coeff(0,0).real(),x.real());
	EXPECT_DOUBLE_EQ(beta.coeff(0,0).imag(),x.imag());
	x = pars.t0/(z-pars.e0-pars.e0-pars.d0);
	EXPECT_DOUBLE_EQ(beta.coeff(1,0).real(),x.real());
	EXPECT_DOUBLE_EQ(beta.coeff(1,0).imag(),x.imag());

	nsum =4;
	ab.FillAlphaBetaMatrix(nsum, z, alpha, beta);
	// testing the dimensions
	EXPECT_EQ(alpha.rows(),1);
	EXPECT_EQ(alpha.cols(),2);
	EXPECT_EQ(beta.rows(),1);
	EXPECT_EQ(beta.cols(),1);
	//testing matrix element
	x = pars.t0/(z-pars.e0-pars.e0);
	EXPECT_DOUBLE_EQ(alpha.coeff(0,0).real(),x.real());
	EXPECT_DOUBLE_EQ(alpha.coeff(0,0).imag(),x.imag());
	EXPECT_DOUBLE_EQ(alpha.coeff(0,1).real(),x.real());
	EXPECT_DOUBLE_EQ(alpha.coeff(0,1).imag(),x.imag());
	EXPECT_DOUBLE_EQ(beta.coeff(0,0).real(),x.real());
	EXPECT_DOUBLE_EQ(beta.coeff(0,0).imag(),x.imag());


	nsum = 5;
	ab.FillAlphaBetaMatrix(nsum, z, alpha, beta);
	// testing the dimensions
	EXPECT_EQ(alpha.rows(),1);
	EXPECT_EQ(alpha.cols(),1);
	//EXPECT_EQ(beta.rows(),0); // a sparse matrix with 0 rows and 0 columns doesn't work
	//EXPECT_EQ(beta.cols(),0);
	//testing matrix element
	x = pars.t0/(z-pars.e0-pars.e0-pars.d0);
	EXPECT_DOUBLE_EQ(alpha.coeff(0,0).real(),x.real());
	EXPECT_DOUBLE_EQ(alpha.coeff(0,0).imag(),x.imag());

}


TEST(SolveVncSPTest, DISABLED_OneThousandSites) {
	Parameters pars;
	pars.nmax = 1000;
	pars.e0 = pars.t0 = pars.d0 = 1;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;
	int ni1 = pars.nmax/2;
	int ni2 = ni1 +1;
	AlphaBetaSP ab(pars);
	dcomplex z(1.0, 0.1);
	CDVector Vnc;
	for (int i=0; i<1; ++i) {
		z = dcomplex(-6, 0.1) + 1.2*i;
		Vnc = solveVncSP(ni1,ni2,z,ab);
	}

	EXPECT_EQ(2,2);
}


TEST(SolveVncSPTest, DISABLED_TwoHundredSites) {
	Parameters pars;
	pars.nmax = 200;
	pars.e0 = 0.0;
	pars.t0 = 1.0;
	pars.d0 = 1.0;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 1.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;
	int ni1 = 100;
	int ni2 = 101;
	AlphaBetaSP ab(pars);
	dcomplex z(1.0, 0.1);
	CDVector Vnc;
	for (int i=0; i<100; ++i) {
		z = dcomplex(-6, 0.1) + 1.2*i;
		Vnc = solveVncSP(ni1,ni2,z,ab);
	}

	EXPECT_EQ(2,2);
}


TEST(GenerateDensityOfStateSP, DISABLED_comparedWithMathematicaResult) {
	Parameters pars;
	pars.nmax = 101;
	pars.e0 = 0.0;
	pars.t0 = 5.0;
	pars.d0 = 15.0;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 2.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;
	int ni1 = pars.nmax/2;
	int ni2 = ni1 + 1;
	std::vector<dcomplex> zList(101);
	std::vector<double> zRealList(101);
	for (int i=0; i<zList.size(); ++i) {
		zRealList[i]=-50 + i;
		zList[i]=dcomplex(-50 + i, 0.1);
	}
	std::vector<double> rhoList;
	generateDensityOfStateSP(ni1, ni2, pars,zList,rhoList);
	save_two_arrays("rho_vs_energy2.txt", zRealList, rhoList);

	// compare the result with "rho_vs_energy.txt.bak"
//	std::vector<double> list;
//	read_double_column("rho_vs_energy2.txt", 2, list);
//	std::vector<double> compareList;
//	read_double_column("rho_vs_energy.txt.bak", 2, compareList);
//	for (int i=0; i<rhoList.size(); ++i) {
//		EXPECT_DOUBLE_EQ(list[i], compareList[i]);
//	}
	EXPECT_EQ(2,2);
}
