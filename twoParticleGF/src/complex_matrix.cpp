/*
 * complexMatrix.hpp
 *
 *  Created on: Dec 19, 2013
 *      Author: pxiang
 */

/*
A simple matrix class (double precision) c++ code

With this class you can:
  - create a 2D matrix with custom size
  - get/set the cell values
  - use operators +, -, *
  - use functions Ones(), Zeros(), Identity(), Inv()
  - print the content of the matrix as an 1D array

Usage:
  you can create a matrix by:
    Matrix A;
    Matrix A = Matrix(rows, cols);
    Matrix A = Matrix(p, rows, cols) where p points to an 1D array of size rows*cols
    Matrix A = B;

  you can get and set matrix elements by:
    A(2,3) = 5.6;    // set an element of Matix A
    value = A(3,1);   // get an element of Matrix A
    value = A.get(3,1); // get an element of a constant Matrix A
    A = B;        // copy content of Matrix B to Matrix A

  you can apply operations with matrices and doubles:
    A = B + C;
    A = B - C;
    A = -B;
    A = B * C;


  the following functions are available:
    A = Ones(rows, cols);
    A = Zeros(rows, cols);
    A = Diag(n);
    A = Inv(B);
    cols = A.GetCols();
    rows = A.GetRows();
    cols = Size(A, 1);
    rows = Size(A, 2);

  you can quick-print the content of a matrix in the console with:
    A.PrintAs1D();
*/


#include "complex_matrix.h"





// calculate the inverse of a general square complex matrix of size N X N
// A points to the array, N is the column size (leading dimension of the array)
void inverse(complex_mkl* A, int N) {
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    complex_mkl *WORK = new complex_mkl[LWORK];
    int INFO;

    zgetrf_(&N,&N,A,&N,IPIV,&INFO);
    zgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    if (INFO!=0) {
    	printf("Failed to find the inverse\n");
    }
    delete IPIV;
    delete WORK;
}


// matrix multiplication A(m,k)*B(k,n) (for complex matrix)
void matrixMul(complex_mkl *A, complex_mkl *B, int m,  int k, int n, complex_mkl *C) {
	complex_mkl alpha = {1.0, 0.0};
	complex_mkl beta = {0.0, 0.0};
    //printf (" Computing matrix product using Intel® MKL dgemm function via CBLAS interface \n\n");
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, &alpha, A, k, B, n, &beta, C, n);
    //printf ("\n Computations completed.\n\n");
}


// matrix multiplication A(m,k)*B(k,n) (real_matrix X complex_matrix)
void matrixMul(double *A, complex_mkl *B, int m, int k, int n,  complex_mkl *C) {
	complex_mkl alpha={1.0, 0.0};
	complex_mkl beta= {0.0, 0.0};
    //printf (" Computing matrix product using Intel® MKL dgemm function via CBLAS interface \n\n");
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, &alpha, A, k, B, n, &beta, C, n);
    //printf ("\n Computations completed.\n\n");
}


// matrix multiplication A(m,k)*B(k,n) (complex_matrix X real_matrix)
void matrixMul(complex_mkl *A, double *B, int m,  int k, int n, complex_mkl *C) {
	complex_mkl alpha = {1.0, 0.0};
	complex_mkl beta = {0.0, 0.0};
    //printf (" Computing matrix product using Intel® MKL dgemm function via CBLAS interface \n\n");
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, &alpha, A, k, B, n, &beta, C, n);
    //printf ("\n Computations completed.\n\n");
}


// solve the linear equation AX = B where A, X, and B are all complex matrices
/*The routine solves the system of linear equations for X:
A*X = B
where
A is a square matrix.
The columns of matrix B are individual right-hand sides.
The columns of X are the corresponding solutions.
The matrix B is overwritten by X. */
void linearSolve(complex_mkl *A, int n, complex_mkl *B, int nrhs) {
	MKL_INT lda = n, ldb = nrhs, info;
	/* Local arrays */
	MKL_INT ipiv[n];
    /* Solve the equations A*X = B */
    info = LAPACKE_zgesv( LAPACK_ROW_MAJOR, n, nrhs, A, lda, ipiv, B, ldb );
    /* Check for the exact singularity */
    if( info > 0 ) {
            printf( "Failed to solve the linear equation AX = B!" );
            exit( 1 );
    }
}


// solve AX = B, the matrix B will be overwritten by X
void solveLinearEqs(ComplexMatrix& A, ComplexMatrix& B) {
	//check the dimensions
	if (A.GetRows()!=B.GetRows()) {
		throw Exception("The matrix dimensions don't match!");
	} else {
		int n = A.GetCols();
		int nrhs = B.GetCols();
		linearSolve(A.GetPointer(), n, B.GetPointer(), nrhs);
	}
}


  // constructor
ComplexMatrix::ComplexMatrix() {
    //printf("Executing constructor Matrix() ...\n");
    // create a Matrix object without content
    p = NULL;
    rows = 0;
    cols = 0;
  }

  // constructor
ComplexMatrix::ComplexMatrix(const int row_count, const int column_count, bool initialize) {
    // create a Matrix object with given number of rows and columns
    p = NULL;

    if (row_count > 0 && column_count > 0){
      rows = row_count;
      cols = column_count;

      p = new complex_mkl[rows*cols];
      if (initialize==true) {
          #pragma omp parallel for
    	  for (int i = 0; i < rows*cols; i++){
    		  // initially fill in zeros for all values in the matrix;
    		  p[i].real = 0.0;
    		  p[i].imag = 0.0;
    	  }
      }
    } else { // if the row_count and column_count is illegal, we just create an empty matrix
    	rows = 0;
    	cols = 0;
    }
  }

  // construct a matrix from an existing one-dimensional array of real numbers
ComplexMatrix::ComplexMatrix(double *a, const int row_count, const int column_count) {
	    // create a Matrix object with given number of rows and columns

	    p = NULL;

	    if (row_count > 0 && column_count > 0){
	      rows = row_count;
	      cols = column_count;

	      p = new complex_mkl[rows*cols];
          #pragma omp parallel for
	      for (int i = 0; i < rows*cols; i++){
	        p[i].real = a[i];
	        p[i].imag = 0.0;
	      }
	    }
  }

  // construct a matrix from an existing one-dimensional array of complex numbers
ComplexMatrix::ComplexMatrix(complex_mkl *a, const int row_count, const int column_count) {
	    // create a Matrix object with given number of rows and columns
	    p = NULL;

	    if (row_count > 0 && column_count > 0){
	      rows = row_count;
	      cols = column_count;

	      p = new complex_mkl[rows*cols];
          #pragma omp parallel for
	      for (int i = 0; i < rows*cols; i++){

	        p[i] = a[i];
	      }
	    }
  }



  // assignment operator
ComplexMatrix::ComplexMatrix(const ComplexMatrix& a){
    rows = a.rows;
    cols = a.cols;
    p = NULL;
    p = new complex_mkl[rows*cols];
    #pragma omp parallel for
    for (int i = 0; i < rows*cols; i++){
      p[i] = a.p[i];
    }
  }

  //

  // index operator. You can use this class like myMatrix(col, row)
  // the indexes are zero based.
  complex_mkl& ComplexMatrix::operator()(const int r, const int c){
    if (p != NULL && r >= 0 && r < rows && c >= 0 && c < cols){
      return p[r*cols + c];
    } else {
      throw Exception("Subscript out of range");
    }
  }

  // index operator. You can use this class like myMatrix.get(col, row)
  // the indexes are zero based.
  // use this function get if you want to read from a const Matrix
  complex_mkl ComplexMatrix::get(const int r, const int c) const {
    if (p != NULL && r >= 0 && r < rows && c >= 0 && c < cols) {
      return p[r*cols + c];
    } else {
      throw Exception("Subscript out of range");
    }
  }


  complex_mkl * ComplexMatrix::GetPointer() {
	  return p;
  }

  complex_mkl * ComplexMatrix::getPointer() const {
	  return p;
  }

  // assignment operator
   ComplexMatrix& ComplexMatrix::operator= (const ComplexMatrix& a) {
     rows = a.rows;
     cols = a.cols;
     if (rows>0 && cols>0) {
    	 delete [] p; //destroy previous memory space
    	 p = NULL;
    	 p = new complex_mkl[rows*cols];
		#pragma omp parallel for
    	 for (int i = 0; i < rows*cols; i++){
    		 p[i] = a.p[i];
    	 }
     } else {
    	 rows = 0;
    	 cols = 0;
    	 if (p!=NULL) {
    		 delete [] p;
    		 p = NULL;
    	 }
     }
     return *this;
   }

  // add a double value (elements wise)
  ComplexMatrix& ComplexMatrix::Add(const double v) {
     #pragma omp parallel for
	 for (int i = 0; i < rows*cols; i++){
		  // initially fill in zeros for all values in the matrix;
		  p[i].real += v;
	  }
	 return *this;
  }

  // add a complex value (elements wise)
  ComplexMatrix& ComplexMatrix::Add(const complex_mkl v) {
     #pragma omp parallel for
	 for (int i = 0; i < rows*cols; i++){
		  // initially fill in zeros for all values in the matrix;
		  p[i].real += v.real;
		  p[i].imag += v.imag;
	  }
	 return *this;
  }


  void ComplexMatrix::AddInPlace(const ComplexMatrix& v) {
	  if (rows!=v.GetRows() && cols!=v.GetCols()) {
		  throw Exception("Dimensions don't match!");
	  }
     #pragma omp parallel for
	 for (int i = 0; i < rows*cols; i++){
		  // initially fill in zeros for all values in the matrix;
		  p[i].real += v.p[i].real;
		  p[i].imag += v.p[i].imag;
	  }
  }

  // subtract a double value (elements wise)
  ComplexMatrix& ComplexMatrix::Subtract(const double v) {
    return Add(-v);
  }

  // subtract a complex value (elements wise)
  ComplexMatrix& ComplexMatrix::Subtract(const complex_mkl v) {
	  complex_mkl v2 = {-v.real, -v.imag};
    return Add(v2);
  }

  // multiply a double value (elements wise)
  ComplexMatrix& ComplexMatrix::Multiply(const double v) {
      #pragma omp parallel for
	  for (int i = 0; i < rows*cols; i++){
		  // initially fill in zeros for all values in the matrix;
		  p[i].real *= v;
	  }
	  return *this;
  }

  // multiply a complex value (elements wise)
  ComplexMatrix& ComplexMatrix::Multiply(const complex_mkl v) {
      #pragma omp parallel for
	  for (int i = 0; i < rows*cols; i++){
		  // initially fill in zeros for all values in the matrix;
		  p[i].real *= v.real;
		  p[i].imag *= v.imag;
	  }
	  return *this;
  }

  // divide a double value (elements wise)
  ComplexMatrix& ComplexMatrix::Divide(const double v) {
     return Multiply(1.0/v);
  }

  // divide a double value (elements wise)
  ComplexMatrix& ComplexMatrix::Divide(const complex_mkl v) {
	 //complex_double w = std::conj(v)/(std::abs(v)*std::abs(v));
	  complex_mkl inverse;
	  inverse.real = v.real/(v.real*v.real + v.imag*v.imag);
	  inverse.imag = -v.imag/(v.real*v.real + v.imag*v.imag);
     return Multiply(inverse);
  }

  // addition of Matrix with Matrix
  ComplexMatrix operator+(const ComplexMatrix& a, const ComplexMatrix& b) {
    // check if the dimensions match
    if (a.rows == b.rows && a.cols == b.cols) {
      ComplexMatrix res(a.rows, a.cols);
	  #pragma omp parallel for
      for (int i = 0; i < a.rows*a.cols; i++) {
         res.p[i].real = a.p[i].real + b.p[i].real;
         res.p[i].imag = a.p[i].imag + b.p[i].imag;
      }
      return res;
    } else {
      // give an error
      throw Exception("Dimensions does not match");
    }

    // return an empty matrix (this never happens but just for safety)
    return ComplexMatrix();
  }




  // addition of Matrix with double
  ComplexMatrix operator+ (const ComplexMatrix& a, const double b) {
    ComplexMatrix res = a;
    res.Add(b);
    return res;
  }

  // addition of Matrix with complex
  ComplexMatrix operator+ (const ComplexMatrix& a, const complex_mkl b) {
    ComplexMatrix res = a;
    res.Add(b);
    return res;
  }

  // addition of double with Matrix
  ComplexMatrix operator+ (const double b, const ComplexMatrix& a) {
    ComplexMatrix res = a;
    res.Add(b);
    return res;
  }


  // subtraction of Matrix with Matrix
  ComplexMatrix operator- (const ComplexMatrix& a, const ComplexMatrix& b) {
    // check if the dimensions match
    if (a.rows == b.rows && a.cols == b.cols) {
      ComplexMatrix res(a.rows, a.cols);
      #pragma omp parallel for
      for (int i = 0; i < a.rows*a.cols; i++) {
         res.p[i].real = a.p[i].real - b.p[i].real;
         res.p[i].imag = a.p[i].imag - b.p[i].imag;
      }
      return res;
    } else {
      // give an error
      throw Exception("Dimensions does not match");
    }

    // return an empty matrix (this never happens but just for safety)
    return ComplexMatrix();
  }




  // subtraction of Matrix with double
  ComplexMatrix operator- (const ComplexMatrix& a, const double b) {
    ComplexMatrix res = a;
    res.Subtract(b);
    return res;
  }

  ComplexMatrix operator- (const ComplexMatrix& a, const complex_mkl b) {
    ComplexMatrix res = a;
    res.Subtract(b);
    return res;
  }

  // subtraction of double with Matrix
  ComplexMatrix operator- (const double b, const ComplexMatrix& a) {
    ComplexMatrix res = -a;
    res.Add(b);
    return res;
  }

  // subtraction of a complex number with Matrix
  ComplexMatrix operator- (const complex_mkl b, const ComplexMatrix& a) {
    ComplexMatrix res = -a;
    res.Add(b);
    return res;
  }

  // operator unary minus
  ComplexMatrix operator- (const ComplexMatrix& a) {
    ComplexMatrix res(a.rows, a.cols);
    #pragma omp parallel for
    for (int i = 0; i < a.rows*a.cols; i++) {
       res.p[i].real = -a.p[i].real;
       res.p[i].imag = -a.p[i].imag;
    }

    return res;
  }

  // operator multiplication
  ComplexMatrix operator* (const ComplexMatrix& a, const ComplexMatrix& b) {
    // check if the dimensions match
    if (a.cols == b.rows) {
      complex_mkl *c = new complex_mkl[a.rows*b.cols];
      //TODO use mkl blas for the matrix multiplication
      matrixMul(a.p, b.p, a.rows, a.cols, b.cols, c);
      ComplexMatrix res = ComplexMatrix(c, a.rows, b.cols);
      delete c;
      return res;
    } else {
      // give an error
      throw Exception("Dimensions does not match");
    }

    // return an empty matrix (this never happens but just for safety)
    return ComplexMatrix();
  }



  // multiplication of Matrix with double
  ComplexMatrix operator* (const ComplexMatrix& a, const double b)
  {
    ComplexMatrix res = a;
    res.Multiply(b);
    return res;
  }
  // multiplication of double with Matrix
  ComplexMatrix operator* (const double b, const ComplexMatrix& a)
  {
    ComplexMatrix res = a;
    res.Multiply(b);
    return res;
  }

  /*
   * returns the inverse of Matrix a
   */
  ComplexMatrix Inv(const ComplexMatrix& a) {
    int rows = a.GetRows();
    int cols = a.GetCols();
    if (rows!=cols) {
  	  throw Exception("Matrix must be square");
    }
    ComplexMatrix res = a;
    inverse(res.p, cols);
    return res;
  }

  // returns the number of rows
  int ComplexMatrix::GetRows() const {
    return rows;
  }

  // returns the number of columns
  int ComplexMatrix::GetCols() const {
    return cols;
  }

  // print the matrix
  void ComplexMatrix::Print() const {
	    printf("\n");
	    for (int r = 0; r < rows; r++) {
	        for (int c = 0; c < cols; c++) {
	            printf("%.5f %+.5f *I  ", p[r*cols + c].real, p[r*cols + c].imag);
	        }
	        printf("\n");
	    }
	    printf("\n");
  }

  // print the content of the matrix as a 1d array
  void ComplexMatrix::PrintAs1D() const {
	  printf("\n");
	  for (int i=0; i<rows*cols; ++i) {
		  printf("(%.5f  %.5f)   ", p[i].real, p[i].imag);
	  }
	  printf("\n");
  }

  // after this call, the matrix will becomes its own inverse
  void ComplexMatrix::InverseInPlace() {
    if (rows!=cols) {
      std::cout << "The matrix is not square!" << std::endl;
  	  throw Exception("The matrix is not square!");
    }

    //inverse(getPointer(), cols);
    inverse(p, cols);
  }


  // destructor
ComplexMatrix::~ComplexMatrix() {
    // clean up allocated memory
    delete [] p;
    p = NULL;
  }


// mixing of real matrix and complex matrix
ComplexMatrix operator+ (const Matrix & a, const ComplexMatrix& b) {
    // check if the dimensions match
	int rows = a.GetRows();
	int cols = a.GetCols();
    if ( rows== b.GetRows() &&  cols== b.GetCols()) {
      ComplexMatrix res(rows, cols);
	  #pragma omp parallel for
      for (int i = 0; i < rows*cols; i++) {
         res.p[i].real = a.getPointer()[i] + b.p[i].real;
         res.p[i].imag = b.p[i].imag;
      }
      return res;
    } else {
      // give an error
      throw Exception("Dimensions does not match");
    }
}

ComplexMatrix operator+ (const ComplexMatrix & a, const Matrix& b) {
    // check if the dimensions match
	int rows = a.GetRows();
	int cols = a.GetCols();
    if ( rows== b.GetRows() &&  cols== b.GetCols()) {
      ComplexMatrix res(rows, cols);
	  #pragma omp parallel for
      for (int i = 0; i < rows*cols; i++) {
         res.p[i].real = a.p[i].real + b.getPointer()[i];
         res.p[i].imag = a.p[i].imag;
      }
      return res;
    } else {
      // give an error
      throw Exception("Dimensions does not match");
    }
}

ComplexMatrix operator- (const Matrix & a, const ComplexMatrix& b) {
    // check if the dimensions match
    if (a.GetRows()== b.GetRows() &&  a.GetCols()== b.GetCols()) {
      int rows = a.GetRows();
      int cols = a.GetCols();
      ComplexMatrix res(rows, cols);
	  #pragma omp parallel for
      for (int r = 0; r < rows; r++) {
    	  for (int c=0; c< cols; c++) {
    		  res(r,c).real = a.get(r,c) - b.get(r,c).real;
    		  res(r,c).imag = -b.get(r,c).imag;
    	  }
      }
      return res;
    } else {
      // give an error
      throw Exception("Dimensions does not match");
    }
}

ComplexMatrix operator- (const ComplexMatrix & a, const Matrix& b) {
    // check if the dimensions match
	int rows = a.GetRows();
	int cols = a.GetCols();
    if ( rows== b.GetRows() &&  cols== b.GetCols()) {
      ComplexMatrix res(rows, cols);
	  #pragma omp parallel for
      for (int i = 0; i < rows*cols; i++) {
         res.p[i].real = a.p[i].real - b.getPointer()[i];
         res.p[i].imag = a.p[i].imag;
      }
      return res;
    } else {
      // give an error
      throw Exception("Dimensions does not match");
    }
}


ComplexMatrix operator* (const Matrix & a, const ComplexMatrix& b) {
	// convert a to a complex matrix
	ComplexMatrix ca = ComplexMatrix(a.getPointer(), a.GetRows(), a.GetCols());
	return operator*(ca, b);
}

ComplexMatrix operator* (const ComplexMatrix & a, const Matrix& b) {
	// convert b to a complex matrix
	ComplexMatrix cb = ComplexMatrix(b.getPointer(), b.GetRows(), b.GetCols());
	return operator*(a, cb);
}



/*******************************************************************************************************/






