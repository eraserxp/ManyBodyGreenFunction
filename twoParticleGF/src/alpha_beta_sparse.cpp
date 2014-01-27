/*
 * alpha_beta_sparse.cpp
 *
 *  Created on: Jan 27, 2014
 *      Author: pxiang
 */

#include "alpha_beta_sparse.h"

// use sparse matrix from eigen c++ library to do the calculations
//
//
//
//
//
//
AlphaBetaSP::AlphaBetaSP(Parameters& ps):ham(ps) {
	nmax = ps.nmax;
	e0 = ps.e0;
	t0 = ps.t0;
	d0 = ps.d0;
	e0MaxDisorder = ps.e0MaxDisorder;
	t0MaxDisorder = ps.t0MaxDisorder;
	d0MaxDisorder = ps.d0MaxDisorder;
	e0seed = ps.e0seed;
	t0seed = ps.t0seed;
	d0seed = ps.d0seed;
	generateIndexMatrix(nmax, Index, VtoG, DimsOfV);
}

void AlphaBetaSP::FillAlphaBetaMatrix(int nsum, dcomplex z, CSMatrix& alpha, CSMatrix& beta) {
	int rows = DimsOfV[nsum];
	int cols = (nsum-1<0)? 0:DimsOfV[nsum-1];
	int cols2 = (nsum+1<2*nmax)? DimsOfV[nsum+1]:0;

	alpha = CSMatrix(rows,cols);
	beta = CSMatrix(rows,cols2);
	TripletList alphaList;
	TripletList betaList;
	alphaList.reserve(max(rows,cols)*2);
	betaList.reserve(max(rows,cols2)*2);

	int i;
	int j;
	dcomplex v_ij, denominator;

	//Eigen::initParallel(); //use multi-threading with eigen c++
	// fill in the matrix row by row
	//#pragma omp parallel for  //-------------------------> may lead to problems
	for (int ith=0; ith<rows; ++ith) {
		Pair p = VtoG(nsum, ith); // obtain g(i,j) which is ith item in V_nsum
		i = p.First();
		j = p.Second();
		denominator = z-ham.energyAtSite(i)-ham.energyAtSite(j)-delta(abs(i-j),1)*ham.d(i, j);
		Neighbors ns = generateNeighbors(i,j, nmax);

		if (ns.firstLeft!=-1) { // the first particle has left neighbor
			int jth = Index(i-1,j);
			//alpha(ith,jth) = ham.t(i-1, i)/denominator;
			v_ij = ham.t(i-1, i)/denominator;
			alphaList.push_back(triplet(ith, jth, v_ij));
		}

		if (ns.firstRight!=-1) { // the first particle has right neighbor
			int jth = Index(i+1,j);
			//beta(ith,jth) = ham.t(i, i+1)/denominator;
			v_ij = ham.t(i, i+1)/denominator;
			betaList.push_back(triplet(ith, jth, v_ij));

		}

		if (ns.secondLeft!=-1) { // the first particle has left neighbor
			int jth = Index(i,j-1);
			//alpha(ith,jth) = ham.t(j-1, j)/denominator;
			v_ij = ham.t(j-1, j)/denominator;
			alphaList.push_back(triplet(ith, jth, v_ij));
		}

		if (ns.secondRight!=-1) { // the first particle has right neighbor
			int jth = Index(i,j+1);
			//beta(ith,jth) = ham.t(j, j+1)/denominator;
			v_ij = ham.t(j, j+1)/denominator;
			betaList.push_back(triplet(ith, jth, v_ij));
		}

	}

	// fill in the matrix
	alpha.setFromTriplets(alphaList.begin(), alphaList.end());
	beta.setFromTriplets(betaList.begin(), betaList.end());
}





int AlphaBetaSP::GetNmax() {
	return nmax;
}

// obtain the factor in front of G(ni1, ni2) in the equations for Green's functions
dcomplex AlphaBetaSP::GetFactor(int ni1, int ni2, dcomplex z) {
	 return z - ham.energyAtSite(ni1)- ham.energyAtSite(ni2)-
				             ham.d(ni1,ni2)*delta(abs(ni1-ni2),1);
}



void changeElements(CSMatrix& csm) {
	int rows = csm.rows();
	int cols = rows;
	for (int k=0; k<csm.outerSize(); ++k) {
	  for (Eigen::SparseMatrix<dcomplex>::InnerIterator it(csm,k); it; ++it) {
	    it.valueRef() = -it.value();
	  }
	}

	for (int i=0; i<rows; ++i) {
		csm.coeffRef(i,i) = 1.0 + csm.coeffRef(i,i);
	}
}


CDMatrix solveDenseLinearEqs(CDMatrix& A, CDMatrix& B) {
	return A.partialPivLu().solve(B);
}


// this requires a lot of conversion between dense matrix and sparse matrix, it is not efficient
CSMatrix solveDenseLinearEqs(CSMatrix& A, CSMatrix& B) {
	CDMatrix AD(A);
	CDMatrix BD(B);
	CDMatrix m= AD.partialPivLu().solve(BD);
	return m.sparseView();
}



CSMatrix solveSparseLinearEqs(CSMatrix& A, CSMatrix& B) {
/******* SparseLU solver ***************/
	A.makeCompressed();
	Eigen::SparseLU< CSMatrix > solver;
/***************************************/

/******* Pardiso solver from intel mkl library *******/
//	Eigen::PardisoLU< CSMatrix > solver;
/****************************************************/

/******* Pardiso solver from intel mkl library *******/
//	Eigen::BiCGSTAB<CSMatrix> solver;
/****************************************************/

	// Compute the ordering permutation vector from the structural pattern of A
	solver.analyzePattern(A);
	// Compute the numerical factorization
	solver.factorize(A);


	CSMatrix X(B.rows(), B.cols());
	//printf("OK 1 \n");
	for (int i=0; i<B.cols(); ++i) {
		CDVector b(B.col(i));
		//printf("OK 2 \n");
		CDVector x = solver.solve(b);
		//printf("OK 3 \n");

/*********** only use matrix elements near the diagonal **********/
//		int range = 5; //10;
//		int start = max(0,i-range);
//		int end = min(i+range+1, x.rows());
//		//printf("OK 4 \n");
//		for (int r=start; r<end; ++r) {
//			X.coeffRef(r,i) = x(r);
//			//printf("OK 5 \n");
//		}
/****************************************************************/

		DVector x_abs = x.array().abs2();
		//double abs_cutoff = x_abs.maxCoeff()/1000.0;
		double abs_cutoff = x_abs.maxCoeff()/1000.0;
		for (int r=0; r<x.rows(); ++r) {
			//if (std::abs(r-i)<=5 || std::abs(r-x.rows())<=5 || x_abs(r) > abs_cutoff) {
			//if (std::abs(r-i)<= 50 || x_abs(r) > abs_cutoff) {
			if (x_abs(r) > abs_cutoff) {

				X.coeffRef(r,i) = x(r);
			}
		}

	}
	//printf("OK 6 \n");
	return X;

//	CSMatrix X = solver.solve(B);
//	if(solver.info()!=Eigen::Success) {
//		printf("Failed to solve AX = B\n");
//	}
//	return result;
}


CDVector solveSparseLinearEqs2(CSMatrix& A, CDVector& b) {
/******* SparseLU solver ***************/
	A.makeCompressed();
	Eigen::SparseLU< CSMatrix > solver;
/***************************************/

/******* Pardiso solver from intel mkl library *******/
//	Eigen::PardisoLU< CSMatrix > solver;
/****************************************************/

/******* Pardiso solver from intel mkl library *******/
//	Eigen::BiCGSTAB<CSMatrix> solver;
/****************************************************/

	// Compute the ordering permutation vector from the structural pattern of A
	solver.analyzePattern(A);
	// Compute the numerical factorization
	solver.factorize(A);

	CDVector x = solver.solve(b);
	if(solver.info()!=Eigen::Success) {
		printf("Failed to solve Ax = b\n");
	}
	return x;
}





CSMatrix fromRightToCenterSP(int Kc, dcomplex z, AlphaBetaSP& ab) {
	int nmax = ab.GetNmax();
	int K = nmax + nmax-1;
	CSMatrix alphan;
	CSMatrix betan;
	ab.FillAlphaBetaMatrix(K,z,alphan, betan);
	CSMatrix AnPlusOne = alphan;
	CSMatrix tmp;

	//printf("Inside fromRightToCenterSP:\n");

	double sparsity;

	for (int n=K-1; n>=Kc+1; n--) {
		//printf("\t nsum = %3d\n", n);
		ab.FillAlphaBetaMatrix(n,z,alphan, betan);
		sparsity = 2.0/alphan.rows();

		tmp = betan*AnPlusOne;
		changeElements(tmp);
		// when the sparsity of alphan and beta is about 1%, switch to the sparse solver
		if (sparsity < 0.01) {
			AnPlusOne = solveSparseLinearEqs(tmp, alphan);
		} else { //use the dense solver
			AnPlusOne = solveDenseLinearEqs(tmp, alphan);
		}
	}

	return AnPlusOne;
}


void checkSparsenessRightToCenter(int Kc, dcomplex z, AlphaBetaSP& ab) {
	int nmax = ab.GetNmax();
	int K = nmax + nmax-1;
	CSMatrix alphan;
	CSMatrix betan;
	ab.FillAlphaBetaMatrix(K,z,alphan, betan);
	CSMatrix AnPlusOne = alphan;
	CSMatrix tmp;
	double sparsity;
	//printf("Inside fromRightToCenterSP:\n");

	for (int n=K-1; n>=Kc+1; n--) {
		ab.FillAlphaBetaMatrix(n,z,alphan, betan);
		sparsity = 2.0/alphan.rows();

		tmp = betan*AnPlusOne;
		changeElements(tmp);
		// when the sparsity of alphan and beta is about 1%, switch to the sparse solver
		if (sparsity < 0.01) {
			AnPlusOne = solveSparseLinearEqs(tmp, alphan);
		} else { //use the dense solver
			AnPlusOne = solveDenseLinearEqs(tmp, alphan);
		}
		double nonZeroPercent = (1.0*AnPlusOne.nonZeros())/(AnPlusOne.rows()*AnPlusOne.cols());
		printf("nsum=%d     non-zero percent:    %5f\n", n, nonZeroPercent);
	}
}



CSMatrix fromLeftToCenterSP(int Kc, dcomplex z, AlphaBetaSP& ab) {
	int nmax = ab.GetNmax();
	int K = 1;
	CSMatrix alphan;
	CSMatrix betan;
	ab.FillAlphaBetaMatrix(K,z,alphan, betan);
	CSMatrix AnMinusOneTilde = betan;
	CSMatrix tmp;
	double sparsity;

	for (int n=K+1; n<=Kc-1; n++) {
		ab.FillAlphaBetaMatrix(n,z,alphan, betan);
		sparsity = 2.0/alphan.rows();

		tmp = alphan*AnMinusOneTilde;
		changeElements(tmp);

		if (sparsity < 0.01) {
			AnMinusOneTilde = solveSparseLinearEqs(tmp, betan);
		} else { //use the dense solver
			AnMinusOneTilde = solveDenseLinearEqs(tmp, betan);
		}
	}

	return AnMinusOneTilde;
}


void checkSparsenessLeftToCenter(int Kc, dcomplex z, AlphaBetaSP& ab) {
	int nmax = ab.GetNmax();
	int K = 1;
	CSMatrix alphan;
	CSMatrix betan;
	ab.FillAlphaBetaMatrix(K,z,alphan, betan);
	CSMatrix AnMinusOneTilde = betan;
	CSMatrix tmp;
	double sparsity;

	for (int n=K+1; n<=Kc-1; n++) {
		ab.FillAlphaBetaMatrix(n,z,alphan, betan);
		sparsity = 2.0/alphan.rows();

		tmp = alphan*AnMinusOneTilde;
		changeElements(tmp);

		if (sparsity < 0.01) {
			AnMinusOneTilde = solveSparseLinearEqs(tmp, betan);
		} else { //use the dense solver
			AnMinusOneTilde = solveDenseLinearEqs(tmp, betan);
		}
		double nonZeroPercent = AnMinusOneTilde.nonZeros();//(100.0*AnMinusOneTilde.nonZeros())/(AnMinusOneTilde.rows()*AnMinusOneTilde.cols());
		nonZeroPercent = nonZeroPercent/(AnMinusOneTilde.rows()*AnMinusOneTilde.cols());
		printf("nsum=%d     row=%d	cols=%d	 non-zero percent:  %5f\n", n, AnMinusOneTilde.rows(),
				                                            AnMinusOneTilde.cols(), nonZeroPercent);
	}

}


void checkSparseness() {
	Parameters pars;
	pars.nmax = 501;
	pars.e0 = 0.0;
	pars.t0 = 5.0;
	pars.d0 = 15.0;
	pars.e0MaxDisorder = pars.t0MaxDisorder = pars.d0MaxDisorder = 2.0;
	pars.e0seed = pars.t0seed = pars.d0seed = 1;
	int ni1 = pars.nmax/2;
	int ni2 = ni1 + 1;
	AlphaBetaSP ab(pars);
	dcomplex z(2.0, 0.1);
	int Kc = ni1 + ni2;
	CSMatrix alphanc;
	CSMatrix betanc;
	ab.FillAlphaBetaMatrix(Kc,z,alphanc, betanc);
	//checkSparsenessRightToCenter(Kc,z,ab);
	checkSparsenessLeftToCenter(Kc,z,ab);
}



CDVector solveVncSP(int ni1, int ni2, dcomplex z, AlphaBetaSP& ab) {
	int Kc = ni1 + ni2;
	CSMatrix alphanc;
	CSMatrix betanc;
	ab.FillAlphaBetaMatrix(Kc,z,alphanc, betanc);
	CSMatrix AncPlusOne = fromRightToCenterSP(Kc,z,ab);
	CSMatrix AncMinusOneTilde = fromLeftToCenterSP(Kc,z,ab);
	CSMatrix mat1 = alphanc*AncMinusOneTilde + betanc*AncPlusOne;
	changeElements(mat1);

	// calculate the dimension of the const vector C
	int nsum = Kc;
	int nth = 0;
	int m;
	for (int i=0; i<= nsum/2; ++i) {
		int j = nsum - i; // i+j = nsum
		if (j>i && j<=ab.GetNmax()) {
			if (i==ni1 && j==ni2) {
				m = nth;
			}
			nth++;
		}
	}

	CDVector constVector = CDVector::Zero(nth);
	dcomplex denominator =  ab.GetFactor(ni1, ni2, z);
	constVector(m) = dcomplex(1.0, 0.0)/denominator;

	return solveSparseLinearEqs2(mat1, constVector);
}






void generateDensityOfStateSP(int ni1, int ni2, Parameters& pars, const std::vector<dcomplex>& zList, std::vector<double>& rhoList) {
	dcomplex z;
	CDVector Vnc;
	AlphaBetaSP ab(pars);
	int nth;
	rhoList.empty();
	//#pragma omp parallel for  ------------> leads to problems
	for (int i=0; i<zList.size(); ++i) {
		z = zList[i];
		Vnc = solveVncSP(ni1,ni2,z,ab);
		nth = getIndex(pars.nmax, ni1+ni2, ni1, ni2);
		rhoList.push_back(-Vnc(nth).imag()/M_PI);
	}
}
