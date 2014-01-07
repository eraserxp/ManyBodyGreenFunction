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
