/*
 * generate_alpha_beta.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: pxiang
 */
#include "alpha_beta.h"

dcomplex convertToDcomplex(complex_mkl& mkl_complex) {
	return dcomplex(mkl_complex.real, mkl_complex.imag);
}

Hamiltonian::Hamiltonian(Parameters& p):rnm_e0(p.nmax+1,1,p.e0seed), rnm_t0(p.nmax+1,p.nmax+1,p.t0seed),
rnm_d0(p.nmax+1,p.nmax+1,p.d0seed){
	nmax = p.nmax;
	e0 = p.e0;
	t0 = p.t0;
	d0 = p.d0;
	e0MaxDisorder = p.e0MaxDisorder;
	t0MaxDisorder = p.t0MaxDisorder;
	d0MaxDisorder = p.d0MaxDisorder;
	e0seed = p.e0seed;
	t0seed = p.t0seed;
	d0seed = p.d0seed;
}

// the energy of excitation at site n when there is no interaction between sites
double Hamiltonian::energyAtSite(int n) {
	return e0 + rnm_e0(n,0)*e0MaxDisorder;
}

// hopping interaction
double Hamiltonian::t(int n, int m) {
	return t0 + rnm_t0(n,m)*t0MaxDisorder;
}

// dynamic interaction
double Hamiltonian::d(int n, int m) {
	return d0 + rnm_d0(n,m)*d0MaxDisorder;
}



// generate an index matrix Index[i][j] = nth which represents the Green's function G(i,j)
// is the nth item in vector V_{i+j}
// we have constricted that G(i,j) for i<j
// DimsOfV: record the dimension of each V_K
// VtoG maps from (i+j, nth) to (i,j)
void generateIndexMatrix(int nmax, IntegerMatrix& Index, PairMatrix& VtoG, std::vector<int>& DimsOfV) {
	int nth;
	int j;
	int nsite = nmax + 1; //total number of sites
	Index = IntegerMatrix(nsite, nsite);
	VtoG = PairMatrix(nmax+nmax, nmax+1);
	DimsOfV.resize(nmax+nmax); //nsum can be as large as nmax+nmax-1
	// set all dimension to be zero initially
	for (int i=0; i<DimsOfV.size(); ++i) {
		DimsOfV[i] = 0;
	}

	for (int nsum = 1; nsum < nmax + nmax; ++nsum) {
		nth = 0;
		for (int i=0; i<= nsum/2; ++i) {
			j = nsum - i; // i+j = nsum
			if (j>i && j<=nmax) {
				Index(i,j) = nth; // g(i,j) is the nth (zero-based) item in V_{i+j} = V_K
				VtoG(nsum,nth) = Pair(i,j); // mapping from (K, nth) to (i,j)
				nth++;
			}
		}
		DimsOfV[nsum] = nth;
	}

}

// given g(ni,nj), it tells the index of g(ni,nj) in vector V_{ni+nj}
// if g(ni,nj) is not in V_{ni+nj}, it returns -1 which means not found
int getIndex(int nmax, int nsum, int ni, int nj) {
	int nth = 0;
	int m;
	if (ni + nj != nsum) {
		return -1;
	}
	for (int i=0; i<= nsum/2; ++i) {
		int j = nsum - i; // i+j = nsum
		if (j>i && j<=nmax) {
			if (i==ni && j==nj) {
				m = nth;
				return m;
			}
			nth++;
		}
	}

	return -1;
}

// obtain the neighbors of the state (i,j)
// we are assuming i<j
Neighbors generateNeighbors(int i, int j, int nmax) {
	 Neighbors ns;
	 if (i-1>=0) {
		 ns.firstLeft = i-1;
	 }

	 if (i+1<=nmax && i+1 != j) {
		 ns.firstRight = i+1;
	 }

	 if (j-1>=0 && j-1 != i) {
		 ns.secondLeft = j-1;
	 }

	 if (j+1<=nmax) {
		 ns.secondRight = j+1;
	 }
	 return ns;
 }


AlphaBeta::AlphaBeta(Parameters& ps):ham(ps) {
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

void AlphaBeta::FillAlphaBetaMatrix(int nsum, complex_mkl z, ComplexMatrix& alpha, ComplexMatrix& beta) {
	int rows = DimsOfV[nsum];
	int cols = (nsum-1<0)? 0:DimsOfV[nsum-1];
	int cols2 = (nsum+1<2*nmax)? DimsOfV[nsum+1]:0;

	alpha = ComplexMatrix(rows,cols,true); // alpha is initialized to 0 initially
	beta = ComplexMatrix(rows,cols2,true);

	// fill in the matrix row by row
	#pragma omp parallel for
	for (int ith=0; ith<rows; ++ith) {
		Pair p = VtoG(nsum, ith); // obtain g(i,j) which is ith item in V_nsum
		int i = p.First();
		int j = p.Second();
		complex_mkl denominator = { z.real -ham.energyAtSite(i)-
				                       ham.energyAtSite(j)-
				                       delta(abs(i-j),1)*ham.d(i, j),
				                       z.imag };
		double denominatorSquared = denominator.real*denominator.real +
				                    denominator.imag*denominator.imag;
		Neighbors ns = generateNeighbors(i,j, nmax);

		if (ns.firstLeft!=-1) { // the first particle has left neighbor
			int jth = Index(i-1,j);
			alpha(ith,jth).real = ham.t(i-1, i)*denominator.real/denominatorSquared;
			alpha(ith,jth).imag = -ham.t(i-1, i)*denominator.imag/denominatorSquared;
		}

		if (ns.firstRight!=-1) { // the first particle has right neighbor
			int jth = Index(i+1,j);
			//beta(ith,jth) = ham.t(i, i+1)/denominator;
			beta(ith,jth).real = ham.t(i, i+1)*denominator.real/denominatorSquared;
			beta(ith,jth).imag = -ham.t(i, i+1)*denominator.imag/denominatorSquared;

		}

		if (ns.secondLeft!=-1) { // the first particle has left neighbor
			int jth = Index(i,j-1);
			//alpha(ith,jth) = ham.t(j-1, j)/denominator;
			alpha(ith,jth).real = ham.t(j-1, j)*denominator.real/denominatorSquared;
			alpha(ith,jth).imag = -ham.t(j-1, j)*denominator.imag/denominatorSquared;
		}

		if (ns.secondRight!=-1) { // the first particle has right neighbor
			int jth = Index(i,j+1);
			//beta(ith,jth) = ham.t(j, j+1)/denominator;
			beta(ith,jth).real = ham.t(j, j+1)*denominator.real/denominatorSquared;
			beta(ith,jth).imag = -ham.t(j, j+1)*denominator.imag/denominatorSquared;
		}

	}
}





int AlphaBeta::GetNmax() {
	return nmax;
}

// obtain the factor in front of G(ni1, ni2) in the equations for Green's functions
complex_mkl AlphaBeta::GetFactor(int ni1, int ni2, complex_mkl z) {
	 return {z.real - ham.energyAtSite(ni1)- ham.energyAtSite(ni2)-
				             ham.d(ni1,ni2)*delta(abs(ni1-ni2),1),
				             z.imag};
}


void changeElements(ComplexMatrix& cm) {
	int rows = cm.GetRows();
	int cols = rows;
	#pragma omp parallel for
	for (int r=0; r<rows; ++r) {
		for (int c=0; c<rows; ++c) {
			cm(r,c).real = -cm(r,c).real;
			cm(r,c).imag = -cm(r,c).imag;
		}
	}

	#pragma omp parallel for
	for (int i=0; i<rows; ++i) {
		cm(i,i).real += 1.0;
	}
}

// The following is a slow version that calculates the inverse of matirx
//ComplexMatrix fromRightToCenter(int Kc, complex_double z, AlphaBeta& ab) {
//	int nmax = ab.GetNmax();
//	int K = nmax + nmax-1;
//	ComplexMatrix alphan;
//	ComplexMatrix betan;
//	ab.FillAlphaBetaMatrix(K,z,alphan, betan);
//	ComplexMatrix AnPlusOne = alphan;
//	ComplexMatrix An;
//	ComplexMatrix tmp1, tmp2;
//
//	for (int n=K-1; n>=Kc+1; n--) {
//		ab.FillAlphaBetaMatrix(n,z,alphan, betan);
//		tmp1 = betan*AnPlusOne;
//		int rows = tmp1.GetRows();
//		tmp2 = Identity(rows) - tmp1;
//		tmp2.InverseInPlace();
//		An = tmp2*alphan;
//		AnPlusOne = An;
//	}
//
//	// try to speed up
//	return An;
//}


ComplexMatrix fromRightToCenter(int Kc, complex_mkl z, AlphaBeta& ab) {
	int nmax = ab.GetNmax();
	int K = nmax + nmax-1;
	ComplexMatrix alphan;
	ComplexMatrix betan;
	ab.FillAlphaBetaMatrix(K,z,alphan, betan);
	ComplexMatrix AnPlusOne = alphan;
	ComplexMatrix tmp1;

	for (int n=K-1; n>=Kc+1; n--) {
		ab.FillAlphaBetaMatrix(n,z,alphan, betan);
		tmp1 = betan*AnPlusOne;
		changeElements(tmp1);
		solveLinearEqs(tmp1, alphan);
		AnPlusOne = alphan;
	}

	return AnPlusOne;
}


// the following is a slow version that calculates the inverse of a matrix explicitly
//ComplexMatrix fromLeftToCenter(int Kc, complex_double z, AlphaBeta& ab) {
//	int nmax = ab.GetNmax();
//	int K = 1;
//	ComplexMatrix alphan;
//	ComplexMatrix betan;
//	ab.FillAlphaBetaMatrix(K,z,alphan, betan);
//	ComplexMatrix AnMinusOneTilde = betan;
//	ComplexMatrix AnTilde;
//	ComplexMatrix tmp1, tmp2;
//
//	for (int n=K+1; n<=Kc-1; n++) {
//		ab.FillAlphaBetaMatrix(n,z,alphan, betan);
//		tmp1 = alphan*AnMinusOneTilde;
//		int rows = tmp1.GetRows();
//		tmp2 = Identity(rows) - tmp1;
//		tmp2.InverseInPlace();
//		AnTilde = tmp2*betan;
//		AnMinusOneTilde = AnTilde;
//	}
//
//	return AnTilde;
//}


ComplexMatrix fromLeftToCenter(int Kc, complex_mkl z, AlphaBeta& ab) {
	int nmax = ab.GetNmax();
	int K = 1;
	ComplexMatrix alphan;
	ComplexMatrix betan;
	ab.FillAlphaBetaMatrix(K,z,alphan, betan);
	ComplexMatrix AnMinusOneTilde = betan;
	ComplexMatrix tmp1;

	for (int n=K+1; n<=Kc-1; n++) {
		ab.FillAlphaBetaMatrix(n,z,alphan, betan);
		tmp1 = alphan*AnMinusOneTilde;
		changeElements(tmp1);
		solveLinearEqs(tmp1, betan); //betan now contains the solution of AX = B
		AnMinusOneTilde = betan;
	}

	return AnMinusOneTilde;
}


//ComplexMatrix solveVnc(int ni1, int ni2, complex_double z, Parameters& pars) {
//	AlphaBeta ab(pars);
ComplexMatrix solveVnc(int ni1, int ni2, complex_mkl z, AlphaBeta& ab) {
	//	AlphaBeta ab(pars);
	// to speed up, we can make ab static --- but it doesn't work
	//static AlphaBeta ab(pars);
	int Kc = ni1 + ni2;
	ComplexMatrix alphanc;
	ComplexMatrix betanc;
	ab.FillAlphaBetaMatrix(Kc,z,alphanc, betanc);
	ComplexMatrix AncPlusOne = fromRightToCenter(Kc,z,ab);
	ComplexMatrix AncMinusOneTilde = fromLeftToCenter(Kc,z,ab);
	ComplexMatrix mat1 = alphanc*AncMinusOneTilde;
	mat1.AddInPlace(betanc*AncPlusOne);
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


	ComplexMatrix constVector(nth, 1, true);
	complex_mkl denominator =  ab.GetFactor(ni1, ni2, z);
	//constVector(m,0) = complex_double(1.0, 0.0)/denominator;
	double square = denominator.real*denominator.real + denominator.imag*denominator.imag;
	constVector(m,0).real = denominator.real/square;
	constVector(m,0).imag = -denominator.imag/square;



	solveLinearEqs(mat1, constVector); // now constVector contains the solution
	return constVector;
}






void generateDensityOfState(int ni1, int ni2, Parameters& pars, const std::vector<complex_mkl>& zList, std::vector<double>& rhoList) {
	complex_mkl z;
	ComplexMatrix Vnc;
	AlphaBeta ab(pars);

	int nth;
	rhoList.empty();
	for (int i=0; i<zList.size(); ++i) {
		z = zList[i];
		Vnc = solveVnc(ni1,ni2,z,ab);
		nth = getIndex(pars.nmax, ni1+ni2, ni1, ni2);
		rhoList.push_back(-Vnc(nth,0).imag/M_PI);
	}
}


















// calculate the density of state for a given real energy
double dos(double E, void* pars) {
	complex_mkl z;
	z.real = E;
	z.imag = 0.1;
	Parameters * parameters = (Parameters *) pars;
	ComplexMatrix Vnc;
	AlphaBeta ab(*parameters);

	int ni1 = parameters->nmax/2;
	int ni2 = ni1 + 1;
	Vnc = solveVnc(ni1,ni2,z,ab);
	int nth = getIndex(parameters->nmax, ni1+ni2, ni1, ni2);
	return -Vnc(nth,0).imag/M_PI;
}

double func_test(double E, void* pars) {

	return sin(E);
}
