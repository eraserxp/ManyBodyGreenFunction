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

#include <list>
#include <vector>


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



class Quartet {

public:
	Quartet(){}

	Quartet(int xi, int yi, int xj, int yj) {
		firstX = xi;
		firstY = yi;
		secondX = xj;
		secondY = yj;
	}

	int FirstX() {
		return firstX;
	}

	int FirstY() {
		return firstY;
	}

	int SecondX() {
		return secondX;
	}

	int SecondY() {
		return secondY;
	}

	//define == for the Quartet class
	bool operator==(const Quartet &other) const {
	  return firstX==other.firstX && firstY==other.firstY
			  && secondX==other.secondX && secondY==other.secondY;
	}


	~Quartet() {}

private:
	int firstX; // x coordinates of the first site
	int firstY;
	int secondX;
	int secondY;
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

typedef std::list<Quartet> QuartetList;
typedef std::vector< Pair > PairVector;

// we use the linked-list to save memory
// each item of the vector is a list of quartet which is expandable
typedef std::vector< QuartetList > QuartetListVector;



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
		firstUp = -1;
		firstDown = -1;
		secondLeft = -1; //the left neighbor of the second particle, -1 means that the left neighbor doesn't exist
		secondRight = -1;
		secondUp = -1;
		secondDown = -1;
	}

	Neighbors2D(int left1, int right1, int up1, int down1, int left2, int right2, int up2, int down2) {
		firstLeft = left1; // the left neighbor of the first particle, -1 means that the left neighbor doesn't exist
		firstRight = right1;
		firstUp = up1;
		firstDown = down1;
		secondLeft = left2; //the left neighbor of the second particle, -1 means that the left neighbor doesn't exist
		secondRight = right2;
		secondUp = up2;
		secondDown = down2;
	}

	bool operator==(const Neighbors2D& other) {
		return firstLeft==other.firstLeft &&
				firstRight==other.firstRight &&
				firstUp==other.firstUp &&
				firstDown==other.firstDown &&
				secondLeft==other.secondLeft &&
				secondRight==other.secondRight &&
				secondUp==other.secondUp &&
				secondDown==other.secondDown;
	}

	int firstLeft;
	int firstRight;
	int firstUp;
	int firstDown;
	int secondLeft;
	int secondRight;
	int secondUp;
	int secondDown;
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

typedef Eigen::MatrixXd DMatrix;

typedef Eigen::ArrayXd DArray;

typedef Eigen::MatrixXi IMatrix; // integer matrix

typedef Eigen::SparseVector< std::complex<double> > CSVector;

typedef Eigen::VectorXcd CDVector;

typedef Eigen::ArrayXcd CDArray;

typedef Eigen::VectorXd DVector;

typedef Eigen::Triplet<std::complex<double> > triplet;
typedef std::vector<triplet> TripletList;


// use boost library
//#include "boost/multi_array.hpp"
//typedef boost::multi_array<int, 4> Array4D;


#endif /* TYPES_H_ */
