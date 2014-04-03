/*
 * matrix_inverse.cpp
 *
 *  Created on: Dec 17, 2013
 *      Author: pxiang
 */

//#include <cstdio>
//
#include "types.h"
#include "mics.h"

#include "random_generator.h"
#include "alpha_beta.h"
#include "alpha_beta_sparse.h"
#include "matrix.h"
#include "complex_matrix.h"
#include "integer_matrix.h"
#include "gtest/gtest.h"
#include "alpha_beta_2D.h"
#include "scattering.h"



void someTest() {
    //using two-dimensional array
    double A2[2][2] = {
         {1,2},
         {3,4}
     };


    complex_mkl B[2*2] = {
    		{1, 0.5}, {2, 1.0},
    		{3,0.2}, {4,0.5}
    };

    inverse(B, 2);
    std::cout << "testing the inverse function (complex matrix)" << std::endl;
    std::cout << "the inverse of the matrix:" << std::endl;
    std::cout << "1.0 + 0.5*I    2.0 + 1.0*I \n3.0 + 0.2*I    4.0 + 0.5*I\n" << "is: " << std::endl;
    printf("%.5f %+.5f*I   %.5f %+.5f*I\n", B[0].real, B[0].imag, B[1].real, B[1].imag);
    printf("%.5f %+.5f*I   %.5f %+.5f*I\n", B[2].real, B[2].imag, B[3].real, B[3].imag);
    printf("\n");

//	  testing calculating the inverse of a very large complex matrix
//    int n = 10000;
//    complex_double *m = new complex_double[n*n];
//    RandomNumberGenerator random;
//
//    for (int i=0; i<n*n; ++i) {
//    	m[i]=random.randomComplex();
//    }
//
//    inverse(m,n);


    // testing the matrix class
    Matrix m = Matrix(2,2, true);
    std::cout << "m = Matrix(2,2, true). m should be all zero:" << std::endl;
    m.Print();

    std::cout<< "After assigning {1,2,3,4} to m, m is now:" << std::endl;
    m(0,0) = 1.0;
    m(0,1) = 2.0;
    m(1,0) = 3.0;
    m(1,1) = 4.0;
    m.Print();

    double A [2*2] = {
        1,2,
        3,4
    };
    Matrix m2 = Matrix(A, 2, 2);
    std::cout << "m2 originally is:" << std::endl;
    m2.Print();


    m2.InverseInPlace();
    std::cout << "After calling InverseInPlace, m2 now is:" << std::endl;
    m2.Print();

    Matrix m3 = Inv(m2);
    std::cout << "Inv(m2) is:" << std::endl;
    m3.Print();

    //testing matrix multiplication
    Matrix m4 = m2*m3;
    std::cout<<"m4 = m2*m3. It should be the identity matrix. Let's print m4:" << std::endl;
    m4.Print();

    // testing the complex matrix class
    ComplexMatrix cm = ComplexMatrix(2,2, true);
    std::cout << "cm should be all zero:" << std::endl;
    cm.Print();

    std::cout<< "After assigning {1,2,3,4} to cm, cm is now:" << std::endl;
    cm(0,0).real = 1.0;
    cm(0,0).imag = 0.0;
    cm(0,1).real = 2.0;
    cm(0,1).imag = 0.0;
    cm(1,0).real = 3.0;
    cm(1,0).imag = 0.0;
    cm(1,1).real = 4.0;
    cm(1,1).imag = 0.0;
    cm.Print();

    ComplexMatrix cm2 = ComplexMatrix(B, 2, 2);
    std::cout << "cm2 originally is:" << std::endl;
    cm2.Print();

    cm2.InverseInPlace();
    std::cout << "After calling InverseInPlace, cm2 now is:" << std::endl;
    cm2.Print();

    ComplexMatrix cm3 = Inv(cm2);
    std::cout << "cm3 = Inv(cm2) is:" << std::endl;
    cm3.Print();

    //testing matrix multiplication
    ComplexMatrix cm4 = cm2*cm3;
    std::cout<<"cm4 = cm2*cm3. It should be the identity matrix. Let's print cm4:" << std::endl;
    cm4.Print();

    //testing matrix multiplication
    cm4 = m2*cm3;
    std::cout << "cm4 = m2*cm3. A real matrix multiplied by a complex matrix. Let's print cm4:" << std::endl;
    cm4.Print();

    cm4 = cm3*m2;
    std::cout << "cm4 = cm3*m2. A complex matrix multiplied by a real matrix. Let's print cm4:" << std::endl;
    cm4.Print();

    // testing integer matrix
    IntegerMatrix im = IntegerMatrix(2,2);
    im = 7;
    std::cout << "An integer matrix im = 7. Let's print im:" << std::endl;
    im.Print();

    //***************************testing PairMatrix*****************************************
    PairMatrix pm(2,2);
    for (int i=0; i<2; ++i) {
    	for (int j=0; j<2; ++j) {
    		pm(i,j) = Pair(2,3);
    	}
    }

    pm(1,0) = Pair(1,0);

    for (int i=0; i<2; ++i) {
    	for (int j=0; j<2; ++j) {
    		printf("pm(%d, %d) = (%d, %d)\n", i, j, pm(i,j).First(), pm(i,j).Second() );
    	}
    }
    //****************************************************************************************

}







void noDynamicInteraction() {
	Parameters pars;
	pars.nmax = 200;
	pars.e0 = 0.0;
	pars.t0 = 5.0;
	pars.d0 = 0.0;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;
	int ni1 = 100;
	int ni2 = 101;
	std::vector<complex_mkl> zList(1001);
	std::vector<double> zRealList(1001);
	for (int i=0; i<zList.size(); ++i) {
		zRealList[i]=-50 + 0.1*i;
		zList[i].real= -50 + 0.1*i;
		zList[i].imag = 0.1;
	}
	std::vector<double> rhoList;
	generateDensityOfState(ni1, ni2, pars,zList,rhoList);
	save_two_arrays("rho_vs_energy_no_dynamic.txt", zRealList, rhoList);
}



void withDisorder() {
	Parameters pars;
	pars.nmax = 200;
	pars.e0 = 0.0;
	pars.t0 = 5.0;
	pars.d0 = 15.0;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 5.0;
	pars.e0seed = 1;
	pars.t0seed = 2;
	pars.d0seed = 3;
	int ni1 = 100;
	int ni2 = 101;
	std::vector<complex_mkl> zList(101);
	std::vector<double> zRealList(101);
	for (int i=0; i<zList.size(); ++i) {
		zRealList[i]=-50 + i;
		zList[i].real= -50 + 0.1*i;
		zList[i].imag = 0.1;
	}
	std::vector<double> rhoList;
	generateDensityOfState(ni1, ni2, pars,zList,rhoList);
	save_two_arrays("rho_vs_energy.txt", zRealList, rhoList);
}


void SpTest() {
	Parameters pars;
	pars.nmax = 200;
	pars.e0 = pars.t0 = pars.d0 = 1;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;

	AlphaBetaSP ab(pars);
	int ni1 = 100;
	int ni2 = 101;
	int Kc = ni1 + ni2;
	dcomplex z(1.0, 0.1);
	CSMatrix AncPlusOne = fromRightToCenterSP(Kc,z,ab);

	CSMatrix ATildeMinusOne = fromLeftToCenterSP(Kc,z,ab);

	CDVector g=solveVncSP(ni1, ni2, z,  ab);
}


void test_2d() {
	Parameters2D pars;
	pars.xmax = 10;
	pars.ymax = 10;
	pars.e0 = 0.0;
	pars.t0 = 1.0;
	pars.d0 = 1.0;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 1.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;
	int x1_i = pars.xmax/2-1;
	int y1_i = pars.ymax/2;
	int x2_i = pars.xmax/2;
	int y2_i = pars.ymax/2;
//	printf("1\n");
	AlphaBeta2D ab(pars);
//	printf("2\n");
	complex_mkl z = {1.0, 0.1};
	ComplexMatrix Vnc;
	for (int i=0; i<100; ++i) {
		z.real = -6 + 1.2*i;
		z.imag = 0.1;
//		printf("3\n");
		Vnc = solveVnc2D(x1_i, y1_i, x2_i, y2_i, z, ab);
//		printf("4\n");
	}

}



void scatteringTest() {
	Parameters pars;
	pars.nmax = 200;
	pars.e0 = 0.0;
	pars.t0 = 1.0;
	pars.d0 = 3.0;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;

	std::vector<double> K_array = linspace(-M_PI*0.99,0.0,10);
//	K_array.push_back(-M_PI/2.0);
//	K_array.push_back(0.0);
//	K_array.push_back(M_PI/2.0);

	std::vector<double> transmission_array;
	int impuritySite = pars.nmax/2;
	double disorderStrength = 15.0;

	for (int i=0; i<K_array.size();++i) {
		ComplexMatrix scatteringState;
		double transmission;
		double K = K_array[i];
		calculateScatteringState(pars, K, impuritySite, disorderStrength, scatteringState, transmission);
		transmission_array.push_back(transmission);
		printf("K=%f     transmission = %f\n", K, transmission);
	}
	// save to the disk
	//save_two_arrays("transmission_vs_k.txt",K_array, transmission_array);

}

void scatteringTest2() {
	Parameters pars;
	pars.nmax = 200;
	pars.e0 = 0.0;
	pars.t0 = 1.0;
	pars.d0 = 5.0;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 0.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;

	std::vector<double> K_array = linspace(-M_PI*0.96,M_PI*0.96,21);
//	K_array.push_back(-M_PI/2.0);
//	K_array.push_back(0.0);
//	K_array.push_back(M_PI/2.0);

	std::vector<double> transmission_array;
	int impuritySite = pars.nmax/2;
	double disorderStrength = 3.5; //13000.0;

	for (int i=0; i<K_array.size();++i) {
		CDVector scatteringState;
		double transmission;
		double K = K_array[i];
		calculateScatteringState2(pars, K, impuritySite, disorderStrength, scatteringState, transmission);
		transmission_array.push_back(transmission);
		printf("K=%f     transmission = %f\n", K, transmission);
	}
	// save to the disk
	save_two_arrays("transmission_vs_k.txt",K_array, transmission_array);

}

//void scattering_direct() {
//	InputParameters pars;
//	pars.nmax = 100;
//	pars.t0 = 1.0;
//	pars.d0 = 5.0;
//	pars.impurityIndex.clear();
//	pars.impurityStrength.clear();
//
//	DVector transmissionCoeff;
//	calculateScatteringState_direct(pars, transmissionCoeff);
//}

int main(int argc, char **argv){

//	Eigen::initParallel(); //use multi-threading with eigen c++

//	someTest();

//	checkSparseness();
//	test_2d();

	scatteringTest2();
//	scatteringTest();

	::testing::InitGoogleTest(&argc, argv);
	int i= RUN_ALL_TESTS();



//	noDynamicInteraction();
	//withDisorder();


//	SpTest();

    return 0;
}
