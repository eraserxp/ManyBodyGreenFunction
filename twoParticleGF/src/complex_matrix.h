/*
 * complex_matrix.h
 *
 *  Created on: Dec 20, 2013
 *      Author: pxiang
 */

#ifndef COMPLEX_MATRIX_H_
#define COMPLEX_MATRIX_H_

#include <cstdlib>
#include <cstdio>
#include <math.h>

#include "types.h"
#include "mkl.h"
#include "mkl_lapacke.h"

#include "matrix.h"

// calculate the inverse of a general square complex matrix of size N X N
// A points to the array, N is the column size (leading dimension of the array)
void inverse(complex_mkl* A, int N);


// matrix multiplication A(m,k)*B(k,n) (for complex matrix)
void matrixMul(complex_mkl *A, complex_mkl *B, int m,  int k, int n, complex_mkl *C);


// matrix multiplication A(m,k)*B(k,n) (real_matrix X complex_matrix)
void matrixMul(double *A, complex_mkl *B, int m, int k, int n,  complex_mkl *C);


// matrix multiplication A(m,k)*B(k,n) (complex_matrix X real_matrix)
void matrixMul(complex_mkl *A, double *B, int m,  int k, int n, complex_mkl *C);

// solve the linear equation AX = B where A, X, and B are all complex matrices
/*The routine solves the system of linear equations for X:
A*X = B
where
A is a square matrix.
The columns of matrix B are individual right-hand sides.
The columns of X are the corresponding solutions.
The matrix B is overwritten by X. */
void linearSolve(complex_mkl *A, int n, complex_mkl *B, int nrhs);



/*
 * Complex Matrix class
 *
 *//*******************************************************************************************************/
class ComplexMatrix {

	  friend class Matrix;

public:
  // constructor
  ComplexMatrix();

  // constructor
  ComplexMatrix(const int row_count, const int column_count, bool initialize=false);

  // construct a matrix from an existing one-dimensional array of real numbers
  ComplexMatrix(double *a, const int row_count, const int column_count);

  // construct a matrix from an existing one-dimensional array of complex numbers
  ComplexMatrix(complex_mkl *a, const int row_count, const int column_count);

  // assignment operator
  //TODO
//  ComplexMatrix(const Matrix& a);


  // assignment operator
  ComplexMatrix(const ComplexMatrix& a);
  //

  // index operator. You can use this class like myMatrix(col, row)
  // the indexes are zero based.
  complex_mkl& operator()(const int r, const int c);

  // index operator. You can use this class like myMatrix.get(col, row)
  // the indexes are zero based.
  // use this function get if you want to read from a const Matrix
  complex_mkl get(const int r, const int c) const;

  complex_mkl * GetPointer();

  complex_mkl * getPointer() const;


  // assignment operator
  //TODO
  //ComplexMatrix& operator= (const Matrix& a);

  // assignment operator
   ComplexMatrix& operator= (const ComplexMatrix& a);

  // add a double value (elements wise)
  ComplexMatrix& Add(const double v);

  // add a complex value (elements wise)
  ComplexMatrix& Add(const complex_mkl v);

  void AddInPlace(const ComplexMatrix& v);


  // subtract a double value (elements wise)
  ComplexMatrix& Subtract(const double v);

  // subtract a complex value (elements wise)
  ComplexMatrix& Subtract(const complex_mkl v);

  // multiply a double value (elements wise)
  ComplexMatrix& Multiply(const double v);

  // multiply a complex value (elements wise)
  ComplexMatrix& Multiply(const complex_mkl v);


  // divide a double value (elements wise)
  ComplexMatrix& Divide(const double v);

  // divide a double value (elements wise)
  ComplexMatrix& Divide(const complex_mkl v);

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
  friend ComplexMatrix operator+(const ComplexMatrix& a, const ComplexMatrix& b);


  // addition of Matrix with double
  friend ComplexMatrix operator+ (const ComplexMatrix& a, const double b) ;

  // addition of Matrix with complex
  friend ComplexMatrix operator+ (const ComplexMatrix& a, const complex_mkl b) ;

  // addition of double with Matrix
  friend ComplexMatrix operator+ (const double b, const ComplexMatrix& a);



  // subtraction of Matrix with Matrix
  friend ComplexMatrix operator- (const ComplexMatrix& a, const ComplexMatrix& b);


  // subtraction of Matrix with double
  friend ComplexMatrix operator- (const ComplexMatrix& a, const double b);

  friend ComplexMatrix operator- (const ComplexMatrix& a, const complex_mkl b);

  // subtraction of a complex number with Matrix
  friend ComplexMatrix operator- (const complex_mkl b, const ComplexMatrix& a);

  // operator unary minus
  friend ComplexMatrix operator- (const ComplexMatrix& a);

  // operator multiplication
  friend ComplexMatrix operator* (const ComplexMatrix& a, const ComplexMatrix& b);


  // multiplication of Matrix with double
  friend ComplexMatrix operator* (const ComplexMatrix& a, const double b);

  // multiplication of double with Matrix
  friend ComplexMatrix operator* (const double b, const ComplexMatrix& a);


  /*
   * returns the inverse of Matrix a
   */
  friend ComplexMatrix Inv(const ComplexMatrix& a);

  // mixing of real matrix and complex matrix
  friend ComplexMatrix operator+ (const Matrix & a, const ComplexMatrix& b);
  friend ComplexMatrix operator+ (const ComplexMatrix & a, const Matrix& b);
  //friend ComplexMatrix operator- (const Matrix & a, const ComplexMatrix& b);
  friend ComplexMatrix operator- (const ComplexMatrix & a, const Matrix& b);
  friend ComplexMatrix operator* (const Matrix & a, const ComplexMatrix& b);
  friend ComplexMatrix operator* (const ComplexMatrix & a, const Matrix& b);




public:
  // destructor
  ~ComplexMatrix();

private:
  int rows;
  int cols;
  complex_mkl * p;     // pointer to a matrix with doubles
};




/*******************************************************************************************************/


ComplexMatrix operator- (const Matrix & a, const ComplexMatrix& b);


// solve AX = B, the matrix B will be overwritten by X
void solveLinearEqs(ComplexMatrix& A, ComplexMatrix& B);


#endif /* COMPLEX_MATRIX_H_ */
