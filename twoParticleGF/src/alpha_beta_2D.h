/*
 * generate_alpha_beta_2D.h
 *
 *  Created on: Jan 6, 2014
 *      Author: pxiang
 */
#ifndef ALPHA_BETA_2D_H_
#define ALPHA_BETA_2D_H_

#include "types.h"
#include "mics.h"
#include "complex_matrix.h"
#include "integer_matrix.h"

#include "random_generator.h"


#include <vector>


// generate the hamiltonian matrix
class Hamiltonian2D {
public:

	Hamiltonian2D(Parameters2D& p);

	double energyAtSite(int nx, int ny);

	double t(int nx, int ny, char direction);

	double d(int nx, int ny, char direction);

private:
	int xmax;
	int ymax;
	double e0, t0, d0;
	double e0MaxDisorder, t0MaxDisorder, d0MaxDisorder;
	unsigned e0seed, t0seed, d0seed;
	OnsiteDisorder  rnm_e0;
	RandomInteraction2D  rnm_t0;
	RandomInteraction2D  rnm_d0;

};


void getCoordinates(int index, int xmax, int ymax, int& xi, int& yi);

int coordinatesToIndex(int xmax, int ymax, int xi, int yi);

#endif
