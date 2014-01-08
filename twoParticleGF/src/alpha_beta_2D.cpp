/*
 * generate_alpha_beta_2D.cpp
 *
 *  Created on: Jan 6, 2014
 *      Author: pxiang
 */

#include "alpha_beta_2D.h"


Hamiltonian2D::Hamiltonian2D(Parameters2D& p):rnm_e0(p.xmax,p.ymax,p.e0seed),
rnm_t0(p.xmax,p.ymax,p.t0seed), rnm_d0(p.xmax,p.ymax,p.d0seed){
	xmax = p.xmax;
	ymax = p.ymax;
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
double Hamiltonian2D::energyAtSite(int nx, int ny) {
	return e0 + rnm_e0(nx,ny)*e0MaxDisorder;
}

// hopping interaction in 4 different direction
double Hamiltonian2D::t(int nx, int ny, char direction) {
	double random = 0.0;
    switch ( direction ) {
       case 'L':
          random = rnm_t0(nx, ny).left;
          break;
       case 'R':
    	  random = rnm_t0(nx, ny).right;
          break;
       case 'U':
    	  random = rnm_t0(nx, ny).up;
          break;
       case 'D':
    	  random = rnm_t0(nx, ny).down;
          break;
       default:
          printf("Illegal direction!");
    }

	return t0 + random*t0MaxDisorder;
}

// dynamic interaction
double Hamiltonian2D::d(int nx, int ny, char direction) {
	double random = 0.0;
    switch ( direction ) {
       case 'L':
          random = rnm_d0(nx, ny).left;
          break;
       case 'R':
    	  random = rnm_d0(nx, ny).right;
          break;
       case 'U':
    	  random = rnm_d0(nx, ny).up;
          break;
       case 'D':
    	  random = rnm_d0(nx, ny).down;
          break;
       default:
          printf("Illegal direction!");
    }

	return t0 + random*d0MaxDisorder;
}

// in 2D, we index the lattice sites in the following way
/*
 *
y

^
|
|
|20	21	22	23	24
|
|15	16	17	18	19
|
|10	11	12	13	14
|
| 5	 6	 7	 8	 9
|
| 0	 1	 2	 3	 4
*---------------------> x
* given the index n for a lattice site i, its coordinates
* are xi = n%cols, yi = n/cols where cols = xmax+1
* for example, the coordinates for site 7 is (7%5, 7/5) = (2, 1)
*/
void getCoordinates(int index, int xmax, int ymax, int& xi, int& yi) {
	xi = index%(xmax + 1);
	yi = index/(xmax + 1);
}

// return -1 if the input coordinates are illegal
int coordinatesToIndex(int xmax, int ymax, int xi, int yi) {
	if (0<=xi && xi<=xmax &&  0<=yi && yi<=ymax) {
		return yi*(xmax + 1) + xi;
	} else {
		return -1;
	}
}


// obtain the neighbors of the state (xi, yi, xj, yj) ---> two sites indexI and indexJ
// to avoid double counting, we only use the lattice state (indexI, indexJ)
// where indexI < indexJ
/*
 * note the neighbors of (indexI, indexJ) should also satisfy the first index < the second index
 */
Neighbors2D generateNeighbors2D(int xi, int yi, int xj, int yj, int xmax, int ymax) {
	 Neighbors2D ns;
	 int indexI = coordinatesToIndex(xmax,ymax, xi, yi);
	 int indexJ = coordinatesToIndex(xmax, ymax, xj, yj);
	 // 4 possible neighbors of i
	 int leftNeighbor = coordinatesToIndex(xmax, ymax, xi-1, yi);
	 if (leftNeighbor!=-1 && leftNeighbor<indexJ) {
		 ns.firstLeft = xi - 1;
	 }

	 int rightNeighbor = coordinatesToIndex(xmax, ymax, xi+1, yi);
	 if (rightNeighbor!=-1 && rightNeighbor<indexJ) {
		 ns.firstRight = xi + 1;
	 }

	 int upNeighbor = coordinatesToIndex(xmax, ymax, xi, yi+1);
	 if (upNeighbor!=-1 && upNeighbor<indexJ) {
		 ns.firstUp = yi + 1;
	 }

	 int downNeighbor = coordinatesToIndex(xmax, ymax, xi, yi-1);
	 if (downNeighbor!=-1 && downNeighbor<indexJ) {
		 ns.firstDown = yi - 1;
	 }



	 // 4 possible neighbors of j
	 leftNeighbor = coordinatesToIndex(xmax, ymax, xj-1, yj);
	 if (leftNeighbor!=-1 && leftNeighbor>indexI) {
		 ns.secondLeft = xj - 1;
	 }

	 rightNeighbor = coordinatesToIndex(xmax, ymax, xj+1, yj);
	 if (rightNeighbor!=-1 && rightNeighbor>indexI) {
		 ns.secondRight = xj + 1;
	 }

	 upNeighbor = coordinatesToIndex(xmax, ymax, xj, yj+1);
	 if (upNeighbor!=-1 && upNeighbor>indexI) {
		 ns.secondUp = yj + 1;
	 }

	 downNeighbor = coordinatesToIndex(xmax, ymax, xj, yj-1);
	 if (downNeighbor!=-1 && downNeighbor>indexI) {
		 ns.secondDown = yj - 1;
	 }

	 return ns;
 }



// generate an index matrix Index[xi][yi][xj][yj] = nth which represents
//the Green's function G(xi, yi, xj, yj) is the nth item in vector V_{xi+yi+xj+yj}
// we have restricted the index for site(xi, yi) must be small than that of (xj, yj)
// DimsOfV: record the dimension of each V_K
// VtoG maps from (xi+yi+xj+yj, nth) to (xi, yi, xj, yj)
void generateIndexMatrix2D(int xmax, int ymax, IntegerMatrix& Index, QuartetListVector& VtoG,
		std::vector<int>& DimsOfV) {
	int nsum;
	//Index = Array4D(boost::extents[xmax+1][ymax+1][xmax+1][ymax+1]);
	int nsite = (xmax+1)*(ymax+1); // total number of sites
	Index = IntegerMatrix(nsite, nsite);
	VtoG.resize(xmax+xmax+ymax+ymax);
	DimsOfV.resize(xmax+xmax+ymax+ymax);


	// set all dimension to be zero initially
	for (int i=0; i<DimsOfV.size(); ++i) {
		DimsOfV[i] = 0;
	}
	// initialize all items of VtoG by an empty list
	for (int i=0; i<DimsOfV.size(); ++i) {
		VtoG[i] = QuartetList();
	}

	for (int x1=0; x1<=xmax; x1++) {
		for (int y1=0; y1<=ymax; y1++) {
			for (int x2=0; x2<=xmax; x2++) {
				for (int y2=0; y2<=ymax; y2++) {
					int index1 = coordinatesToIndex(xmax,ymax,x1,y1);
					int index2 = coordinatesToIndex(xmax,ymax,x2,y2);
					if (index1 < index2) {
						nsum = x1 + y1 + x2 + y2;
						Index(index1, index2) = DimsOfV[nsum];
						DimsOfV[nsum]+=1;
						VtoG[nsum].push_back(Quartet(x1,y1,x2,y2));
					}
				}
			}
		}
	}

}




AlphaBeta2D::AlphaBeta2D(Parameters2D& ps):ham(ps) {
	xmax = ps.xmax;
	ymax = ps.ymax;
	e0 = ps.e0;
	t0 = ps.t0;
	d0 = ps.d0;
	e0MaxDisorder = ps.e0MaxDisorder;
	t0MaxDisorder = ps.t0MaxDisorder;
	d0MaxDisorder = ps.d0MaxDisorder;
	e0seed = ps.e0seed;
	t0seed = ps.t0seed;
	d0seed = ps.d0seed;
	generateIndexMatrix2D(xmax, ymax, Index, VtoG, DimsOfV);
}

void AlphaBeta2D::FillAlphaBetaMatrix(int nsum, complex_mkl z, ComplexMatrix& alpha, ComplexMatrix& beta) {
	int rows = DimsOfV[nsum];
	int cols = (nsum-1>=1)? DimsOfV[nsum-1]:0; //nsum-1 can't smaller than 1
	int cols2 = (nsum+1<=xmax+ymax+xmax+ymax-1)? DimsOfV[nsum+1]:0; // nsum+1 can't be larger than that
	// if cols and cols exceed the range, we just set them to zero
//	printf("1\n");
	alpha = ComplexMatrix(rows,cols,true); // alpha is initialized to 0 initially
//	printf("2\n");
	beta = ComplexMatrix(rows,cols2,true);
//	printf("3\n");

	// fill in the matrix row by row
	#pragma omp parallel for
	for (QuartetList::iterator it=VtoG[nsum].begin(); it!=VtoG[nsum].end(); ++it) {
//		printf("4\n");
		Quartet p = *it; // retrieve the coordinates (xi,yi,xj,yj) from the vector V_nsum
		int x1 = p.FirstX(); // coordinates for the first particle
		int y1 = p.FirstY();
		int x2 = p.SecondX(); // coordinates for the second particle
		int y2 = p.SecondY();
//		printf("5\n");
		int index1 = coordinatesToIndex(xmax,ymax,x1,y1); // index for the first particle
		int index2 = coordinatesToIndex(xmax,ymax,x2,y2); // index for the second particle
		int ith = Index(index1,index2); // the order of state (index1, index2) in V_nsum
		complex_mkl denominator = { z.real -ham.energyAtSite(x1, y1)-
				                       ham.energyAtSite(x2, y2)-
				                   // first particle is the left nearest neighbor of the second
				                   delta(y1-y2,0)*delta(x2-x1,1)*ham.d(x2, y2, 'L') -
				                   // first particle is the down nearest neighbor of the second
				                   delta(x1-x2,0)*delta(y2-y1,1)*ham.d(x2, y2, 'D'),
				                       z.imag };
		double denominatorSquared = denominator.real*denominator.real +
				                    denominator.imag*denominator.imag;
		// generate all the neighbors for the two sites index1 and index2
		Neighbors2D ns = generateNeighbors2D(x1, y1, x2, y2, xmax, ymax);

		// consider the four possible neighbors of particle 1
		if (ns.firstLeft!=-1) { // the first particle has left neighbor
			// index for the left neighbor of particle 1
			int index1_left = coordinatesToIndex(xmax,ymax,x1-1,y1);
			int jth = Index(index1_left,index2);
			double t = ham.t(x1, y1, 'L');
			alpha(ith,jth).real = t*denominator.real/denominatorSquared;
			alpha(ith,jth).imag = -t*denominator.imag/denominatorSquared;
		}

		if (ns.firstRight!=-1) { // the first particle has right neighbor
			// index for the right neighbor of particle 1
			int index1_right = coordinatesToIndex(xmax,ymax,x1+1,y1);
			int jth = Index(index1_right,index2);
			double t = ham.t(x1, y1, 'R');
			beta(ith,jth).real = t*denominator.real/denominatorSquared;
			beta(ith,jth).imag = -t*denominator.imag/denominatorSquared;

		}

		if (ns.firstUp!=-1) { // the first particle has up neighbor
			// index for the up neighbor of particle 1
			int index1_up = coordinatesToIndex(xmax,ymax,x1,y1+1);
			int jth = Index(index1_up,index2);
			double t = ham.t(x1, y1, 'U');
			beta(ith,jth).real = t*denominator.real/denominatorSquared;
			beta(ith,jth).imag = -t*denominator.imag/denominatorSquared;
		}

		if (ns.firstDown!=-1) { // the first particle has down neighbor
			// index for the up neighbor of particle 1
			int index1_down = coordinatesToIndex(xmax,ymax,x1,y1-1);
			int jth = Index(index1_down,index2);
			double t = ham.t(x1, y1, 'D');
			alpha(ith,jth).real = t*denominator.real/denominatorSquared;
			alpha(ith,jth).imag = -t*denominator.imag/denominatorSquared;
		}


		// consider the four possible neighbors of particle 2
		if (ns.secondLeft!=-1) { // the second particle has left neighbor
			// index for the left neighbor of particle 2
			int index2_left = coordinatesToIndex(xmax,ymax,x2-1,y2);
			int jth = Index(index1,index2_left);
			double t = ham.t(x2, y2, 'L');
			alpha(ith,jth).real = t*denominator.real/denominatorSquared;
			alpha(ith,jth).imag = -t*denominator.imag/denominatorSquared;
		}

		if (ns.secondRight!=-1) { // the second particle has right neighbor
			// index for the right neighbor of particle 2
			int index2_right = coordinatesToIndex(xmax,ymax,x2+1,y2);
			int jth = Index(index1,index2_right);
			double t = ham.t(x2, y2, 'R');
			beta(ith,jth).real = t*denominator.real/denominatorSquared;
			beta(ith,jth).imag = -t*denominator.imag/denominatorSquared;

		}

		if (ns.secondUp!=-1) { // the second particle has up neighbor
			// index for the up neighbor of particle 2
			int index2_up = coordinatesToIndex(xmax,ymax,x2,y2+1);
			int jth = Index(index1,index2_up);
			double t = ham.t(x2, y2, 'U');
			beta(ith,jth).real = t*denominator.real/denominatorSquared;
			beta(ith,jth).imag = -t*denominator.imag/denominatorSquared;
		}

		if (ns.secondDown!=-1) { // the second particle has down neighbor
			// index for the up neighbor of particle 2
			int index2_down = coordinatesToIndex(xmax,ymax,x2,y2-1);
			int jth = Index(index1,index2_down);
			double t = ham.t(x2, y2, 'D');
			alpha(ith,jth).real = t*denominator.real/denominatorSquared;
			alpha(ith,jth).imag = -t*denominator.imag/denominatorSquared;
		}


	}
}


int AlphaBeta2D::GetIndex(int x1, int y1, int x2, int y2) {
	int index1 = coordinatesToIndex(xmax,ymax,x1,y1);
	int index2 = coordinatesToIndex(xmax,ymax,x2,y2);
	return Index(index1, index2);
}

int AlphaBeta2D::GetDimOfV(int nsum) {
	return DimsOfV[nsum];
}

int AlphaBeta2D::GetXmax() {
	return xmax;
}

int AlphaBeta2D::GetYmax() {
	return ymax;
}

// obtain the factor in front of G(x1i, y1i, x2i, y2i) in the equations for Green's functions
complex_mkl AlphaBeta2D::GetFactor(int x1_i, int y1_i, int x2_i, int y2_i, complex_mkl z) {
	 return {z.real - ham.energyAtSite(x1_i, y1_i)- ham.energyAtSite(x2_i, y2_i)-
		 // first particle is the left nearest neighbor of the second
		 	delta(y1_i - y2_i, 0)*delta(x2_i - x1_i, 1)*ham.d(x2_i, y2_i, 'L') -
		 // first particle is the down nearest neighbor of the second
		 	delta(x1_i - x2_i, 0)*delta(y2_i - y1_i, 1)*ham.d(x2_i, y2_i, 'D'),
				             z.imag};
}


ComplexMatrix fromRightToCenter2D(int Kc, complex_mkl z, AlphaBeta2D& ab) {
	int xmax = ab.GetXmax();
	int ymax = ab.GetYmax();
	int K = xmax + ymax + xmax + ymax -1;
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


ComplexMatrix fromLeftToCenter2D(int Kc, complex_mkl z, AlphaBeta2D& ab) {
	int xmax = ab.GetXmax();
	int ymax = ab.GetYmax();
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
ComplexMatrix solveVnc2D(int x1_i, int y1_i, int x2_i, int y2_i, complex_mkl z, AlphaBeta2D& ab) {
	int Kc = x1_i + y1_i + x2_i + y2_i;
	ComplexMatrix alphanc;
	ComplexMatrix betanc;
//	printf("1\n");
	ab.FillAlphaBetaMatrix(Kc,z,alphanc, betanc);
//	printf("2\n");
	ComplexMatrix AncPlusOne = fromRightToCenter2D(Kc,z,ab);
//	printf("3\n");
	ComplexMatrix AncMinusOneTilde = fromLeftToCenter2D(Kc,z,ab);
//	printf("4\n");
	ComplexMatrix mat1 = alphanc*AncMinusOneTilde;
	mat1.AddInPlace(betanc*AncPlusOne);
	changeElements(mat1);

	// calculate the dimension of the const vector C
	int n = ab.GetDimOfV(Kc);
	// find out which term is not zero
	int m = ab.GetIndex(x1_i, y1_i, x2_i, y2_i);

	ComplexMatrix constVector(n, 1, true);
	complex_mkl denominator =  ab.GetFactor(x1_i, y1_i, x2_i, y2_i, z);
	//constVector(m,0) = complex_double(1.0, 0.0)/denominator;
	double square = denominator.real*denominator.real + denominator.imag*denominator.imag;
	constVector(m,0).real = denominator.real/square;
	constVector(m,0).imag = -denominator.imag/square;

	solveLinearEqs(mat1, constVector); // now constVector contains the solution
	return constVector;
}


void generateDensityOfState2D(int x1_i, int y1_i, int x2_i, int y2_i, Parameters2D& pars,
		        const std::vector<complex_mkl>& zList, std::vector<double>& rhoList) {
	complex_mkl z;
	ComplexMatrix Vnc;
	AlphaBeta2D ab(pars);

	int nth;
	rhoList.empty();
	for (int i=0; i<zList.size(); ++i) {
		z = zList[i];
		Vnc = solveVnc2D(x1_i, y1_i, x2_i, y2_i,z,ab);
		nth = ab.GetIndex(x1_i, y1_i, x2_i, y2_i);
		rhoList.push_back(-Vnc(nth,0).imag/M_PI);
	}
}
