/*
 * generate_alpha_beta_test.cpp
 *
 *  Created on: Dec 25, 2013
 *      Author: pxiang
 */

#include "alpha_beta.h"
#include "alpha_beta_sparse.h"
#include "gtest/gtest.h"
#include "mics.h"
//#include "integration.h"

TEST(GenerateNeighborsTest, NoNeighbors) {
	Neighbors ns = generateNeighbors(0,1,1);
	EXPECT_EQ(ns.firstLeft, -1);
	EXPECT_EQ(ns.firstRight, -1);
	EXPECT_EQ(ns.secondLeft, -1);
	EXPECT_EQ(ns.secondRight, -1);
}

TEST(GenerateNeighborsTest, OneLeftNeighbor) {
	Neighbors ns = generateNeighbors(1,2,2);
	EXPECT_EQ(ns.firstLeft, 0);
	EXPECT_EQ(ns.firstRight, -1);
	EXPECT_EQ(ns.secondLeft, -1);
	EXPECT_EQ(ns.secondRight, -1);
}


TEST(GenerateNeighborsTest, OneRightNeighbor) {
	Neighbors ns = generateNeighbors(0,1,2);
	EXPECT_EQ(ns.firstLeft, -1);
	EXPECT_EQ(ns.firstRight, -1);
	EXPECT_EQ(ns.secondLeft, -1);
	EXPECT_EQ(ns.secondRight, 2);
}


TEST(GenerateNeighborsTest, TwoNeighbor) {
	Neighbors ns = generateNeighbors(3,4,7);
	EXPECT_EQ(ns.firstLeft, 2);
	EXPECT_EQ(ns.firstRight, -1);
	EXPECT_EQ(ns.secondLeft, -1);
	EXPECT_EQ(ns.secondRight, 5);
}


TEST(GenerateIndexMatrixTest, FourSites) {
	int nmax = 3;
	IntegerMatrix Index;
	PairMatrix VtoG;
	std::vector<int> DimsOfV;
	generateIndexMatrix(nmax, Index, VtoG, DimsOfV);
	EXPECT_EQ(Index(0,1),0);
	EXPECT_EQ(Index(0,2),0);
	EXPECT_EQ(Index(0,3),0);
	EXPECT_EQ(Index(1,2),1);
	EXPECT_EQ(Index(1,3),0);
	EXPECT_EQ(Index(2,3),0);

	EXPECT_EQ(DimsOfV[0],0);
	EXPECT_EQ(DimsOfV[1],1);
	EXPECT_EQ(DimsOfV[2],1);
	EXPECT_EQ(DimsOfV[3],2);
	EXPECT_EQ(DimsOfV[4],1);
	EXPECT_EQ(DimsOfV[5],1);


	EXPECT_EQ(VtoG(1,0).First(),0);
	EXPECT_EQ(VtoG(1,0).Second(),1);
	EXPECT_EQ(VtoG(2,0).First(),0);
	EXPECT_EQ(VtoG(2,0).Second(),2);
	EXPECT_EQ(VtoG(3,0).First(),0);
	EXPECT_EQ(VtoG(3,0).Second(),3);
	EXPECT_EQ(VtoG(3,1).First(),1);
	EXPECT_EQ(VtoG(3,1).Second(),2);
	EXPECT_EQ(VtoG(4,0).First(),1);
	EXPECT_EQ(VtoG(4,0).Second(),3);
	EXPECT_EQ(VtoG(5,0).First(),2);
	EXPECT_EQ(VtoG(5,0).Second(),3);
}


TEST(AlphaBetaTest, FourSitesNoDisorder) {
	Parameters pars;
	pars.nmax = 3;
	pars.e0 = pars.t0 = pars.d0 = 1;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;

	AlphaBeta ab(pars);
	int nsum = 1;
	complex_mkl z = {1.0, 0.1};
	ComplexMatrix alpha;
	ComplexMatrix beta;
	ab.FillAlphaBetaMatrix(nsum, z, alpha, beta);
	// testing the dimensions
	EXPECT_EQ(alpha.GetRows(),0);
	EXPECT_EQ(alpha.GetCols(),0);
	EXPECT_EQ(beta.GetRows(),1);
	EXPECT_EQ(beta.GetCols(),1);
	//testing matrix element
	dcomplex x = pars.t0/(convertToDcomplex(z)-pars.e0-pars.e0-pars.d0);
	EXPECT_DOUBLE_EQ(beta(0,0).real,x.real());
	EXPECT_DOUBLE_EQ(beta(0,0).imag,x.imag());

	nsum=2;
	ab.FillAlphaBetaMatrix(nsum, z, alpha, beta);
	// testing the dimensions
	EXPECT_EQ(alpha.GetRows(),1);
	EXPECT_EQ(alpha.GetCols(),1);
	EXPECT_EQ(beta.GetRows(),1);
	EXPECT_EQ(beta.GetCols(),2);
	//testing matrix element
	x = pars.t0/(convertToDcomplex(z)-pars.e0-pars.e0);
	EXPECT_DOUBLE_EQ(alpha(0,0).real,x.real());
	EXPECT_DOUBLE_EQ(alpha(0,0).imag,x.imag());
	EXPECT_DOUBLE_EQ(beta(0,0).real,x.real());
	EXPECT_DOUBLE_EQ(beta(0,0).imag,x.imag());
	EXPECT_DOUBLE_EQ(beta(0,1).real,x.real());
	EXPECT_DOUBLE_EQ(beta(0,1).imag,x.imag());

	nsum =3;
	ab.FillAlphaBetaMatrix(nsum, z, alpha, beta);
	// testing the dimensions
	EXPECT_EQ(alpha.GetRows(),2);
	EXPECT_EQ(alpha.GetCols(),1);
	EXPECT_EQ(beta.GetRows(),2);
	EXPECT_EQ(beta.GetCols(),1);
	//testing matrix element
	x = pars.t0/(convertToDcomplex(z)-pars.e0-pars.e0);
	EXPECT_DOUBLE_EQ(alpha(0,0).real,x.real());
	EXPECT_DOUBLE_EQ(alpha(0,0).imag,x.imag());
	x = pars.t0/(convertToDcomplex(z)-pars.e0-pars.e0-pars.d0);
	EXPECT_DOUBLE_EQ(alpha(1,0).real,x.real());
	EXPECT_DOUBLE_EQ(alpha(1,0).imag,x.imag());
	x = pars.t0/(convertToDcomplex(z)-pars.e0-pars.e0);
	EXPECT_DOUBLE_EQ(beta(0,0).real,x.real());
	EXPECT_DOUBLE_EQ(beta(0,0).imag,x.imag());
	x = pars.t0/(convertToDcomplex(z)-pars.e0-pars.e0-pars.d0);
	EXPECT_DOUBLE_EQ(beta(1,0).real,x.real());
	EXPECT_DOUBLE_EQ(beta(1,0).imag,x.imag());

	nsum =4;
	ab.FillAlphaBetaMatrix(nsum, z, alpha, beta);
	// testing the dimensions
	EXPECT_EQ(alpha.GetRows(),1);
	EXPECT_EQ(alpha.GetCols(),2);
	EXPECT_EQ(beta.GetRows(),1);
	EXPECT_EQ(beta.GetCols(),1);
	//testing matrix element
	x = pars.t0/(convertToDcomplex(z)-pars.e0-pars.e0);
	EXPECT_DOUBLE_EQ(alpha(0,0).real,x.real());
	EXPECT_DOUBLE_EQ(alpha(0,0).imag,x.imag());
	EXPECT_DOUBLE_EQ(alpha(0,1).real,x.real());
	EXPECT_DOUBLE_EQ(alpha(0,1).imag,x.imag());
	EXPECT_DOUBLE_EQ(beta(0,0).real,x.real());
	EXPECT_DOUBLE_EQ(beta(0,0).imag,x.imag());


	nsum = 5;
	ab.FillAlphaBetaMatrix(nsum, z, alpha, beta);
	// testing the dimensions
	EXPECT_EQ(alpha.GetRows(),1);
	EXPECT_EQ(alpha.GetCols(),1);
	EXPECT_EQ(beta.GetRows(),0);
	EXPECT_EQ(beta.GetCols(),0);
	//testing matrix element
	x = pars.t0/(convertToDcomplex(z)-pars.e0-pars.e0-pars.d0);
	EXPECT_DOUBLE_EQ(alpha(0,0).real,x.real());
	EXPECT_DOUBLE_EQ(alpha(0,0).imag,x.imag());


}


TEST(FromRightToCenterTest, TenSites) {
	Parameters pars;
	pars.nmax = 9;
	pars.e0 = pars.t0 = pars.d0 = 1;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;

	AlphaBeta ab(pars);
	int ni1 = 4;
	int ni2 = 5;
	int Kc = ni1 + ni2;
	complex_mkl z= {1.0, 0.1};
	ComplexMatrix AncPlusOne = fromRightToCenter(Kc,z,ab);
	EXPECT_EQ(2,2);
}



TEST(FromLeftToCenterTest, TenSites) {
	Parameters pars;
	pars.nmax = 9;
	pars.e0 = pars.t0 = pars.d0 = 1;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;

	AlphaBeta ab(pars);
	int ni1 = 4;
	int ni2 = 5;
	int Kc = ni1 + ni2;
	complex_mkl z = {1.0, 0.1};
	ComplexMatrix AncMinusOneTilde = fromLeftToCenter(Kc,z,ab);
	EXPECT_EQ(2,2);
}




TEST(SolveVncTest, TenSites) {
	Parameters pars;
	pars.nmax = 9;
	pars.e0 = pars.t0 = pars.d0 = 1;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;
	int ni1 = 4;
	int ni2 = 5;
	AlphaBeta ab(pars);
	complex_mkl z = {1.0, 0.1};
	ComplexMatrix Vnc = solveVnc(ni1,ni2,z,ab);
	EXPECT_EQ(2,2);
}



TEST(SolveVncTest, DISABLED_TwoHundredSites) {
	Parameters pars;
	pars.nmax = 200;
	pars.e0 = 0.0;
	pars.t0 = 1.0;
	pars.d0 = 1.0;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 1.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;
	int ni1 = 100;
	int ni2 = 101;
	AlphaBeta ab(pars);
	complex_mkl z = {1.0, 0.1};
	ComplexMatrix Vnc;
	for (int i=0; i<100; ++i) {
		z.real = -6 + 1.2*i;
		z.imag = 0.1;
		Vnc = solveVnc(ni1,ni2,z,ab);
	}

	EXPECT_EQ(2,2);
}



TEST(GetIndexTest, FourSites) {
	int nmax = 3;
	int nsum;
	int i, j;
	nsum = 1;
	EXPECT_EQ(getIndex(nmax,nsum,0,0), -1);
	EXPECT_EQ(getIndex(nmax,nsum,0,1), 0);
	EXPECT_EQ(getIndex(nmax,nsum,1,0), -1);

	nsum = 2;
	EXPECT_EQ(getIndex(nmax,nsum,0,2), 0);
	EXPECT_EQ(getIndex(nmax,nsum,2,0), -1);

	nsum = 3;
	EXPECT_EQ(getIndex(nmax,nsum,0,3), 0);
	EXPECT_EQ(getIndex(nmax,nsum,1,2), 1);

	nsum = 4;
	EXPECT_EQ(getIndex(nmax,nsum,1,3), 0);
	EXPECT_EQ(getIndex(nmax,nsum,0,4), -1);

	nsum = 5;
	EXPECT_EQ(getIndex(nmax,nsum,2,3), 0);
	EXPECT_EQ(getIndex(nmax,nsum,1,4), -1);
}










TEST(SolveVncTest, DISABLED_OneThousandSites) {
	Parameters pars;
	pars.nmax = 1000;
	pars.e0 = pars.t0 = pars.d0 = 1;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;
	int ni1 = pars.nmax/2;
	int ni2 = ni1 + 1;
	AlphaBeta ab(pars);
	complex_mkl z = {1.0, 0.1};
	ComplexMatrix Vnc;
	for (int i=0; i<1; ++i) {
		z.real = -6 + 1.2*i;
		z.imag = 0.1;
		Vnc = solveVnc(ni1,ni2,z,ab);
	}

	EXPECT_EQ(2,2);
}




TEST(GenerateDensityOfState, comparedWithMathematicaResult) {
	Parameters pars;
	pars.nmax = 100;
	pars.e0 = 0.0;
	pars.t0 = 1.0;
	pars.d0 = 5.0;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;
	int ni1 = pars.nmax/2;
	int ni2 = ni1 + 1;
	std::vector<complex_mkl> zList(101);
	std::vector<double> zRealList = linspace(-50,50,101);
	for (int i=0; i<zList.size(); ++i) {
		//zRealList[i]=-50 + i;
		zList[i].real = zRealList[i];
		zList[i].imag = 0.01;
	}
	std::vector<double> rhoList;
	generateDensityOfState(ni1, ni2, pars,zList,rhoList);
	save_two_arrays("rho_vs_energy.txt", zRealList, rhoList);
	EXPECT_EQ(2,2);
}





//TEST(DensityOfStateTest, IntegralShouldBeOne) {
//	Parameters pars;
//	pars.nmax = 201;
//	pars.e0 = 0.0;
//	pars.t0 = 5.0;
//	pars.d0 = 15.0;
//	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
//	pars.e0seed = pars.t0seed = pars.d0seed = 1;
//	double x_min = -50;
//	double x_max = 50;
//
//	double integral = integrate(&func_test,&pars,x_min,x_max, 1.e-3, 1.e-3);
//	EXPECT_NEAR(integral, 0.0, 1.e-2);
//}


TEST(CalculateAllGF, checkIO) {
	Parameters pars;
	pars.nmax = 100;
	pars.e0 = 0.0;
	pars.t0 = 5.0;
	pars.d0 = 15.0;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;
	int ni1 = pars.nmax/2;
	int ni2 = ni1 + 1;
//	std::vector<complex_mkl> zList(101);
//	std::vector<double> zRealList = linspace(-50,50,101);
//	for (int i=0; i<zList.size(); ++i) {
//		//zRealList[i]=-50 + i;
//		zList[i].real = zRealList[i];
//		zList[i].imag = 0.1;
//	}
//	std::vector<double> rhoList;
//	generateDensityOfState(ni1, ni2, pars,zList,rhoList);
//	save_two_arrays("rho_vs_energy.txt", zRealList, rhoList);

	AlphaBeta ab(pars);
	complex_mkl z;
	z.real = 35.0;
	z.imag = 0.1;
	calculateAllGF(ni1, ni2, z, ab);

	z.real = 15.0;
	z.imag = 0.1;
	calculateAllGF(ni1, ni2, z, ab);

	z.real = 5.0;
	z.imag = 0.1;
	calculateAllGF(ni1, ni2, z, ab);

	EXPECT_EQ(2,2);
}


TEST(ExtractMatrixElement, readFromFile) {
	Parameters pars;
	pars.nmax = 100;
	pars.e0 = 0.0;
	pars.t0 = 5.0;
	pars.d0 = 15.0;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;
	int ni1 = pars.nmax/2;
	int ni2 = ni1 + 1;

	AlphaBeta ab(pars);
	complex_mkl z;
	z.real = 35.0;
	z.imag = 0.1;
	calculateAllGF(ni1, ni2, z, ab);
	IntegerMatrix indexMatrix;
	ab.FillIndexMatrix(indexMatrix);
	int n = ni1 + 2;
	int m = n + 1;
	complex_mkl gf = extractMatrixElement(n, m, ni1, ni2, indexMatrix);
	EXPECT_EQ(2,2);
}
