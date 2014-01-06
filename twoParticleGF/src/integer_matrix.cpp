/*
 * integer_matrix.cpp
 *
 *  Created on: Dec 23, 2013
 *      Author: pxiang
 */

  // constructor

#include "integer_matrix.h"



IntegerMatrix::IntegerMatrix() {
    //printf("Executing constructor Matrix() ...\n");
    // create a Matrix object without content
    p = NULL;
    rows = 0;
    cols = 0;
  }

  // constructor
IntegerMatrix::IntegerMatrix(const int row_count, const int column_count, bool initialize) {
    // create a Matrix object with given number of rows and columns
    p = NULL;

    if (row_count > 0 && column_count > 0){
      rows = row_count;
      cols = column_count;

      p = new int[rows*cols];
      if (initialize==true) {
          #pragma omp parallel for
    	  for (int i = 0; i < rows*cols; i++){
    		  // initially fill in zeros for all values in the matrix;
    		  p[i] = 0;
    	  }
      }
    }
  }

  // construct a matrix from an existing one-dimensional array of real numbers
IntegerMatrix::IntegerMatrix(int *a, const int row_count, const int column_count) {
	    // create a Matrix object with given number of rows and columns
	    p = NULL;

	    if (row_count > 0 && column_count > 0){
	      rows = row_count;
	      cols = column_count;

	      p = new int[rows*cols];
          #pragma omp parallel for
	      for (int i = 0; i < rows*cols; i++){
	        p[i] = a[i];
	      }
	    }
  }

// index operator. You can use this class like myMatrix(col, row)
 // the indexes are zero based.
 int& IntegerMatrix::operator()(const int r, const int c){
   if (p != NULL && r >= 0 && r < rows && c >= 0 && c < cols){
     return p[r*cols + c];
   } else {
     throw Exception("Subscript out of range");
   }
 }


 // assignment operator
  IntegerMatrix& IntegerMatrix::operator= (const IntegerMatrix& a) {
    rows = a.rows;
    cols = a.cols;
    p = new int[rows*cols];
    #pragma omp parallel for
    for (int i = 0; i < rows*cols; i++){
      p[i] = a.p[i];
    }
    return *this;
  }


  // assignment operator, use an integer to initialize the Integer Matrix
  // before calling it, the space of matrix must be allocated
   IntegerMatrix& IntegerMatrix::operator= (int a) {
     #pragma omp parallel for
     for (int i = 0; i < rows*cols; i++){
       p[i] = a;
     }
     return *this;
   }

   // print the matrix
   void IntegerMatrix::Print() const {
 	    printf("\n");
 	    for (int r = 0; r < rows; r++) {
 	        for (int c = 0; c < cols; c++) {
 	            printf("%d  ", p[r*cols + c]);
 	        }
 	        printf("\n");
 	    }
 	    printf("\n");
   }

	// print the content of the matrix as a 1d array
	void IntegerMatrix::PrintAs1D() const {
		  printf("\n");
		  for (int i=0; i<rows*cols; ++i) {
			  printf("%d  ", p[i]);
		  }
		  printf("\n");
	}


   // destructor
   IntegerMatrix::~IntegerMatrix() {
     // clean up allocated memory
     delete [] p;
     p = NULL;
   }
