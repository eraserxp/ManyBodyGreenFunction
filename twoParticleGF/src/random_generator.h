/*
 * random_generator.h
 *
 *  Created on: Dec 17, 2013
 *      Author: pxiang
 */

#ifndef RANDOM_GENERATOR_H_
#define RANDOM_GENERATOR_H_

#include "types.h"

class RandomNumberGenerator {
	private:
		unsigned seed;

	public:


// constructor
	RandomNumberGenerator(unsigned inputSeed) {
		seed = inputSeed;
		srand(seed);
	}

	// default constructor
	RandomNumberGenerator() {
		//use the current time as seed for random numbers
		seed = (unsigned)time(NULL);
		srand(seed);
	}

	void SetSeed(unsigned inputSeed) {
		seed = inputSeed;
		srand(seed);
	}

// produce a complex random number in [0, 1) + [0, 1)*I given the seed
	complex_mkl randomComplex() {
		double randomReal = (double)rand()/(double)RAND_MAX;
		double randomImag = (double)rand()/(double)RAND_MAX;
		complex_mkl rc = {randomReal,randomImag};
		return rc;
	}

	// produce a complex random number in [0, 1) + [0, 1)*I given the seed
	double randomReal() {
		double random = (double)rand()/(double)RAND_MAX;
		return random;
	}


};

// generate a matrix of given dimension whose elements are random numbers in [0, 1)
class RandomNumberMatrix {
public:
	RandomNumberMatrix(int row_count, int col_count, unsigned seed) {
		rows = row_count;
		cols = col_count;
		p = NULL;
		p = new double[rows*cols];
		rng.SetSeed(seed);
		for (int i=0; i<rows*cols; i++) {
			p[i] = rng.randomReal();
		}
	}

	double operator()(const int r, const int c) const {
		return p[r*cols + c];
	}


	~RandomNumberMatrix() {
		delete [] p;
		p = NULL;
		rows = 0;
		cols = 0;
	}

private:
	double *p;
	int rows;
	int cols;
	RandomNumberGenerator rng;
};

class OnsiteDisorder {
public:
	OnsiteDisorder(int xmax, int ymax, unsigned seed) {
		rows = xmax;
		cols = ymax;
		p = new double[rows*cols];
		rng.SetSeed(seed);
		for (int i=0; i<rows*cols; i++) {
			p[i] = rng.randomReal();
		}
	}

	double operator()(const int nx, const int ny) const {
		return p[nx*cols + ny];
	}

private:
	int rows;
	int cols;
	double *p;
	RandomNumberGenerator rng;
};


// generate a matrix of given dimension whose elements are random numbers in [0, 1)
class RandomInteraction2D {
public:
	RandomInteraction2D(int row_count, int col_count, unsigned seed) {
		rows = row_count;
		cols = col_count;
		p = NULL;
		p = new InteractionInFourDirections[rows*cols];
		rng.SetSeed(seed);
		for (int i=0; i<rows*cols; i++) {
			p[i].left = rng.randomReal();
			p[i].right = rng.randomReal();
			p[i].up = rng.randomReal();
			p[i].down = rng.randomReal();
		}
	}

	InteractionInFourDirections operator()(const int r, const int c) const {
		return p[r*cols + c];
	}


	~RandomInteraction2D() {
		delete [] p;
		p = NULL;
		rows = 0;
		cols = 0;
	}

private:
	InteractionInFourDirections *p;
	int rows;
	int cols;
	RandomNumberGenerator rng;
};

#endif
