/*
 * types.h
 *
 *  Created on: Dec 17, 2013
 *      Author: pxiang
 */

#ifndef TYPES_H_
#define TYPES_H_

#include <math.h>
#define _USE_MATH_DEFINES

#include<iostream>

#include <complex>

#include "mkl.h"

typedef std::complex<double> dcomplex;

//user-defined mkl complex type
//#define MKL_Complex16 complex_double
typedef MKL_Complex16 complex_mkl;



// define my own pair
class Pair {

public:
	Pair(){}

	Pair(int i, int j) {
		first = i;
		second = j;
	}

	int First() {
		return first;
	}

	int Second() {
		return second;
	}

	~Pair() {}

private:
	int first;
	int second;
};

// two-dimensional array of Pairs
class PairMatrix {
public:
	// constructor
	PairMatrix() {
		p = NULL;
		rows = 0;
		cols = 0;
	}

	// constructor
	PairMatrix(int rows_count, int cols_count) {
		rows = rows_count;
		cols = cols_count;
		p = new Pair*[rows];
		for (int i=0; i<rows; ++i) {
			p[i] = new Pair[cols];
		}
	}

	//assignment operator
	PairMatrix& operator= (const PairMatrix& pm) {
		// destroy the original one
		if (p!=NULL) {
			for (int i=0; i<rows; ++i) {
				delete [] p[i];
				p[i]= NULL;
			}
			delete [] p;
			p = NULL;
		}
		// construct a new one
		rows = pm.rows;
		cols = pm.cols;
		p = new Pair*[rows];
		for (int i=0; i<rows; ++i) {
			p[i] = new Pair[cols];
		}
		// assign the value of pm to p
		for (int i=0; i<rows; ++i) {
			for (int j=0; j<cols; ++j) {
				p[i][j] = pm.p[i][j];
			}
		}
		return *this;
	}

	// return the matrix element m(i,j)
	Pair& operator()(int i, int j) {
		return p[i][j];
	}

	// destructor
	~PairMatrix() {
		if (p!=NULL) {
			for (int i=0; i<rows; ++i) {
				delete [] p[i];
				p[i]= NULL;
			}
			delete [] p;
			p = NULL;
		}
	}

private:
	int rows;
	int cols;
	Pair **p;
};


class Neighbors {
public:
	Neighbors() {
		firstLeft = -1; // the left neighbor of the first particle, -1 means that the left neighbor doesn't exist
		firstRight = -1;
		secondLeft = -1; //the left neighbor of the second particle, -1 means that the left neighbor doesn't exist
		secondRight = -1;
	}

	int firstLeft;
	int firstRight;
	int secondLeft;
	int secondRight;
};


// describe the possible 8 neighbors in 2D (for two particles)
class Neighbors2D {
public:
	Neighbors2D() {
		firstLeft = -1; // the left neighbor of the first particle, -1 means that the left neighbor doesn't exist
		firstRight = -1;
		firstTop = -1;
		firstBottom = -1;
		secondLeft = -1; //the left neighbor of the second particle, -1 means that the left neighbor doesn't exist
		secondRight = -1;
		secondTop = -1;
		secondBottom = -1;
	}

	int firstLeft;
	int firstRight;
	int firstTop;
	int firstBottom;
	int secondLeft;
	int secondRight;
	int secondTop;
	int secondBottom;
};


// used to describe the interaction strengths between a particle
// and its 4 neighours in 4 directions
class InteractionInFourDirections {
public:
	InteractionInFourDirections() {
		left = 0.0;
		right = 0.0;
		up = 0.0;
		down = 0.0;
	}

	double left;
	double right;
	double up;
	double down;
};


typedef struct {
	int nmax;
	double e0, t0, d0;
	double e0MaxDisorder, t0MaxDisorder, d0MaxDisorder;
	unsigned e0seed, t0seed, d0seed;
} Parameters;


typedef struct {
	int xmax; // maximum index in x axis (starting from zero)
	int ymax;// maximum index in y axis (starting from zero)
	double e0, t0, d0;
	double e0MaxDisorder, t0MaxDisorder, d0MaxDisorder;
	unsigned e0seed, t0seed, d0seed;
} Parameters2D;


// use eigen c++ library
//#define EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport> // use intel sparse solver
#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<std::complex<double> > CSMatrix;

typedef Eigen::MatrixXcd CDMatrix;

typedef Eigen::SparseVector< std::complex<double> > CSVector;

typedef Eigen::VectorXcd CDVector;

typedef Eigen::VectorXd DVector;

typedef Eigen::Triplet<std::complex<double> > triplet;
typedef std::vector<triplet> TripletList;



#endif /* TYPES_H_ */
