/*
 * integer_matrix.h
 *
 *  Created on: Dec 23, 2013
 *      Author: pxiang
 */

#ifndef INTEGER_MATRIX_H_
#define INTEGER_MATRIX_H_

#include <cstdlib>
#include <cstdio>
#include <math.h>
#include "matrix.h"

class IntegerMatrix {


public:
  // constructor
	IntegerMatrix();

  // constructor
	IntegerMatrix(const int row_count, const int column_count, bool initialize=false);

  // construct a matrix from an existing one-dimensional array of integer numbers
	IntegerMatrix(int *a, const int row_count, const int column_count);

	int& operator()(const int r, const int c);

	IntegerMatrix& operator= (const IntegerMatrix& a);

	IntegerMatrix& operator= (int a);

	void Print() const;

	// print the content of the matrix as a 1d array
	void PrintAs1D() const;

  // destructor
  ~IntegerMatrix();

private:
  int rows;
  int cols;
  int * p;     // pointer to a matrix with doubles
};

#endif /* INTEGER_MATRIX_H_ */
