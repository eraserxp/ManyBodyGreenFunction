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

// hopping interaction
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

	 int topNeighbor = coordinatesToIndex(xmax, ymax, xi, yi+1);
	 if (topNeighbor!=-1 && topNeighbor<indexJ) {
		 ns.firstTop = yi + 1;
	 }

	 int bottomNeighbor = coordinatesToIndex(xmax, ymax, xi, yi-1);
	 if (bottomNeighbor!=-1 && bottomNeighbor<indexJ) {
		 ns.firstBottom = yi - 1;
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

	 topNeighbor = coordinatesToIndex(xmax, ymax, xj, yj+1);
	 if (topNeighbor!=-1 && topNeighbor>indexI) {
		 ns.secondTop = yj + 1;
	 }

	 bottomNeighbor = coordinatesToIndex(xmax, ymax, xj, yj-1);
	 if (bottomNeighbor!=-1 && bottomNeighbor>indexI) {
		 ns.secondBottom = yj - 1;
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
