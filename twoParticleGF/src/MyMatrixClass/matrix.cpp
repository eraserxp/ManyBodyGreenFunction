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




#include "matrix.h"


// calculate the inverse of a general double matrix
void inverse(double* A, int N) {
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    if (INFO!=0) {
    	printf("Failed to find the inverse\n");
    }

    delete IPIV;
    delete WORK;
}

// matrix multiplication A(m,k)*B(k,n)
void matrixMul(double *A, double *B, int m, int k, int n,  double *C) {
	double alpha = 1.0;
	double beta = 0.0;
    //printf (" Computing matrix product using Intel® MKL dgemm function via CBLAS interface \n\n");
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, alpha, A, k, B, n, beta, C, n);
    //printf ("\n Computations completed.\n\n");
}










  Matrix::Matrix() {
    //printf("Executing constructor Matrix() ...\n");
    // create a Matrix object without content
    p = NULL;
    rows = 0;
    cols = 0;
  }

  // constructor
  Matrix::Matrix(const int row_count, const int column_count, bool initialize) {
    // create a Matrix object with given number of rows and columns
    p = NULL;

    if (row_count > 0 && column_count > 0){
      rows = row_count;
      cols = column_count;

      p = new double[rows*cols];
      if (initialize==true) {
          #pragma omp parallel for
    	  for (int i = 0; i < rows*cols; i++){
    		  // initially fill in zeros for all values in the matrix;
    		  p[i] = 0;
    	  }
      }
    }
  }

  // construct a matrix from an existing one-dimensional array
  Matrix::Matrix(double *a, const int row_count, const int column_count) {
	    // create a Matrix object with given number of rows and columns
	    p = NULL;

	    if (row_count > 0 && column_count > 0){
	      rows = row_count;
	      cols = column_count;

	      p = new double[rows*cols];
          #pragma omp parallel for
	      for (int i = 0; i < rows*cols; i++){

	        p[i] = a[i];
	      }
	    }
  }

  // assignment operator
  Matrix::Matrix(const Matrix& a){
    rows = a.rows;
    cols = a.cols;
    p = new double[rows*cols];
    #pragma omp parallel for
    for (int i = 0; i < rows*cols; i++){
      // initially fill in zeros for all values in the matrix;
      p[i] = a.p[i];
    }
  }

  //

  // index operator. You can use this class like myMatrix(col, row)
  // the indexes are zero based.
  double& Matrix::operator()(const int r, const int c){
    if (p != NULL && r >= 0 && r < rows && c >= 0 && c < cols){
      return p[r*cols + c];
    } else {
      throw Exception("Subscript out of range");
    }
  }

  // index operator. You can use this class like myMatrix.get(col, row)
  // the indexes are zero based.
  // use this function get if you want to read from a const Matrix
  double Matrix::get(const int r, const int c) const {
    if (p != NULL && r >= 0 && r < rows && c >= 0 && c < cols) {
      return p[r*cols + c];
    } else {
      throw Exception("Subscript out of range");
    }
  }

  double * Matrix::getPointer() const {
	  return p;
  }

  // get the element as if the matrix is an 1d array
  double Matrix::get1D(const int i) const {
	  return p[i];
  }


  // assignment operator
  Matrix& Matrix::operator= (const Matrix& a) {
    rows = a.rows;
    cols = a.cols;
    p = new double[rows*cols];
    #pragma omp parallel for
    for (int i = 0; i < rows*cols; i++){
      // initially fill in zeros for all values in the matrix;
      p[i] = a.p[i];
    }
    return *this;
  }

  // add a double value (elements wise)
  Matrix& Matrix::Add(const double v) {
     #pragma omp parallel for
	 for (int i = 0; i < rows*cols; i++){
		  // initially fill in zeros for all values in the matrix;
		  p[i] += v;
	  }
	 return *this;
  }

  // subtract a double value (elements wise)
  Matrix& Matrix::Subtract(const double v) {
    return Add(-v);
  }

  // multiply a double value (elements wise)
  Matrix& Matrix::Multiply(const double v) {
      #pragma omp parallel for
	  for (int i = 0; i < rows*cols; i++){
		  // initially fill in zeros for all values in the matrix;
		  p[i] *= v;
	  }
	  return *this;
  }

  // divide a double value (elements wise)
  Matrix& Matrix::Divide(const double v) {
     return Multiply(1/v);
  }


  // returns the number of rows
  int Matrix::GetRows() const {
    return rows;
  }

  // returns the number of columns
  int Matrix::GetCols() const {
    return cols;
  }

  // print the matrix
  void Matrix::Print() const {
	    printf("\n");
	    for (int r = 0; r < rows; r++) {
	        for (int c = 0; c < cols; c++) {
	            printf("%.5f  ", p[r*cols + c]);
	        }
	        printf("\n");
	    }
	    printf("\n");
  }


  // print the content of the matrix as a 1d array
  void Matrix::PrintAs1D() const {
	  printf("\n");
	  for (int i=0; i<rows*cols; ++i) {
		  printf("%.5f  ", p[i]);
	  }
	  printf("\n");
  }

  // after this call, the matrix will becomes its own inverse
  void Matrix::InverseInPlace() {
    if (rows!=cols) {
      std::cout << "The matrix is not square!" << std::endl;
  	  throw Exception("The matrix is not square!");
    }

    //inverse(getPointer(), cols);
    inverse(p, cols);
  }

  Matrix::~Matrix() {
	  delete [] p;
	  p = NULL;
  }


  // addition of Matrix with Matrix
  Matrix operator+(const Matrix& a, const Matrix& b) {
    // check if the dimensions match
    if (a.rows == b.rows && a.cols == b.cols) {
      Matrix res(a.rows, a.cols);
	  #pragma omp parallel for
      for (int i = 0; i < a.rows*a.cols; i++) {
         res.p[i] = a.p[i] + b.p[i];
      }
      return res;
    } else {
      // give an error
      throw Exception("Dimensions does not match");
    }

    // return an empty matrix (this never happens but just for safety)
    return Matrix();
  }

  // addition of Matrix with double
  Matrix operator+ (const Matrix& a, const double b) {
    Matrix res = a;
    res.Add(b);
    return res;
  }

  // addition of double with Matrix
  Matrix operator+ (const double b, const Matrix& a) {
    Matrix res = a;
    res.Add(b);
    return res;
  }

  // subtraction of Matrix with Matrix
  Matrix operator- (const Matrix& a, const Matrix& b) {
    // check if the dimensions match
    if (a.rows == b.rows && a.cols == b.cols) {
      Matrix res(a.rows, a.cols);
      #pragma omp parallel for
      for (int i = 0; i < a.rows*a.cols; i++) {
         res.p[i] = a.p[i] - b.p[i];
      }
      return res;
    } else {
      // give an error
      throw Exception("Dimensions does not match");
    }

    // return an empty matrix (this never happens but just for safety)
    return Matrix();
  }

  // subtraction of Matrix with double
  Matrix operator- (const Matrix& a, const double b) {
    Matrix res = a;
    res.Subtract(b);
    return res;
  }

  // subtraction of double with Matrix
  Matrix operator- (const double b, const Matrix& a) {
    Matrix res = -a;
    res.Add(b);
    return res;
  }

  // operator unary minus
  Matrix operator- (const Matrix& a) {
    Matrix res(a.rows, a.cols);
    #pragma omp parallel for
    for (int i = 0; i < a.rows*a.cols; i++) {
       res.p[i] = -a.p[i];
    }

    return res;
  }

  // operator multiplication
  Matrix operator* (const Matrix& a, const Matrix& b) {
    // check if the dimensions match
    if (a.cols == b.rows) {
      double *c = new double[a.rows*b.cols];
      //TODO use mkl blas for the matrix multiplication
      matrixMul(a.p, b.p, a.rows, a.cols, b.cols, c);
      Matrix res = Matrix(c, a.rows, b.cols);
      delete c;
      return res;
    } else {
      // give an error
      throw Exception("Dimensions does not match");
    }

    // return an empty matrix (this never happens but just for safety)
    return Matrix();
  }

  // multiplication of Matrix with double
  Matrix operator* (const Matrix& a, const double b) {
    Matrix res = a;
    res.Multiply(b);
    return res;
  }
  // multiplication of double with Matrix
  Matrix operator* (const double b, const Matrix& a) {
    Matrix res = a;
    res.Multiply(b);
    return res;
  }

  /*
   * returns the inverse of Matrix a
   */
  Matrix Inv(const Matrix& a) {
    int rows = a.GetRows();
    int cols = a.GetCols();
    if (rows!=cols) {
  	  throw Exception("Matrix must be square");
    }
    Matrix res = a;
    inverse(res.p, cols);
    return res;
  }



/**
 * returns a matrix with size cols x rows with ones as values
 */
Matrix Ones(const int rows, const int cols) {
  Matrix res = Matrix(rows, cols);
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      res(r, c) = 1;
    }
  }
  return res;
}

/**
 * returns a matrix with size cols x rows with zeros as values
 */
Matrix Zeros(const int rows, const int cols) {
  return Matrix(rows, cols);
}


/**
 * returns an identity matrix with size n x n
 */
Matrix Identity(const int n) {
  Matrix res = Matrix(n, n, true);
  #pragma omp parallel for
  for (int i = 0; i < n; i++) {
    res(i, i) = 1;
  }
  return res;
}





