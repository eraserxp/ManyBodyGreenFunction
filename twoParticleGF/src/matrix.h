/*
 * matrix.h
 *
 *  Created on: Dec 20, 2013
 *      Author: pxiang
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <iostream>

//#include "types.h"
#include "mkl.h"


// Declarations
void inverse(double* A, int N);
void matrixMul(double *A, double *B, int m, int k, int n,  double *C);

/*
 * a simple exception class
 * you can create an exeption by entering
 *   throw Exception("...Error description...");
 * and get the error message from the data msg for displaying:
 *   err.msg
 */
class Exception {
public:
  const char* msg;
  Exception(const char* arg)
   : msg(arg)
  {
  }
};

class Matrix {



public:
  // constructor
  Matrix();

  // constructor
  Matrix(const int row_count, const int column_count, bool initialize=false);

  // construct a matrix from an existing one-dimensional array
  Matrix(double *a, const int row_count, const int column_count);

  // assignment operator
  Matrix(const Matrix& a);

  //

  // index operator. You can use this class like myMatrix(col, row)
  // the indexes are zero based.
  double& operator()(const int r, const int c);

  // index operator. You can use this class like myMatrix.get(col, row)
  // the indexes are zero based.
  // use this function get if you want to read from a const Matrix
  double get(const int r, const int c) const ;

  double * getPointer() const;

  // get the element as if the matrix is an 1d array
  double get1D(const int i) const ;

  // assignment operator
  Matrix& operator= (const Matrix& a);

  // add a double value (elements wise)
  Matrix& Add(const double v);

  // subtract a double value (elements wise)
  Matrix& Subtract(const double v);

  // multiply a double value (elements wise)
  Matrix& Multiply(const double v);

  // divide a double value (elements wise)
  Matrix& Divide(const double v);

  // returns the number of rows
  int GetRows() const;

  // returns the number of columns
  int GetCols() const;

  void Print() const;

  // print the content of the matrix as a 1d array
  void PrintAs1D() const;

  // after this call, the matrix will becomes its own inverse
  void InverseInPlace();

  // addition of Matrix with Matrix
  friend Matrix operator+(const Matrix& a, const Matrix& b);

  // addition of Matrix with double
  friend Matrix operator+ (const Matrix& a, const double b);
  // addition of double with Matrix
  friend Matrix operator+ (const double b, const Matrix& a);

  // subtraction of Matrix with Matrix
  friend Matrix operator- (const Matrix& a, const Matrix& b);

  // subtraction of Matrix with double
  friend Matrix operator- (const Matrix& a, const double b);

  // subtraction of double with Matrix
  friend Matrix operator- (const double b, const Matrix& a);

  // operator unary minus
  friend Matrix operator- (const Matrix& a);

  // operator multiplication
  friend Matrix operator* (const Matrix& a, const Matrix& b);

  // multiplication of Matrix with double
  friend Matrix operator* (const Matrix& a, const double b);

  // multiplication of double with Matrix
  friend Matrix operator* (const double b, const Matrix& a);

  /*
   * returns the inverse of Matrix a
   */
  friend Matrix Inv(const Matrix& a);


public:
  // destructor
  ~Matrix();

private:
  int rows;
  int cols;
  double * p;     // pointer to a matrix with doubles
};


Matrix Identity(const int n);
Matrix Inv(const Matrix& a);
Matrix Ones(const int rows, const int cols);
Matrix Zeros(const int rows, const int cols);


#endif /* MATRIX_H_ */
