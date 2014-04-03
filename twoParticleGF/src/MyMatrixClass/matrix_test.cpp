/*
 * matrix_test.cpp
 *
 *  Created on: Dec 24, 2013
 *      Author: pxiang
 */

#include "matrix.h"
#include "gtest/gtest.h"

TEST(InverseTest, Real2X2Matrix) {
    double A [2*2] = {
        1,2,
        3,4
    };


    inverse(A, 2);
    EXPECT_DOUBLE_EQ(A[0], -2.0) << "The calculated inverse is wrong!";
    EXPECT_DOUBLE_EQ(A[1], 1.0) << "The calculated inverse is wrong!";
    EXPECT_DOUBLE_EQ(A[2], 1.5) << "The calculated inverse is wrong!";
    EXPECT_DOUBLE_EQ(A[3], -0.5) << "The calculated inverse is wrong!";

}


TEST(MatrixMulTest, RealMatrix) {
    double A [2*2] = {
        1,2,
        3,4
    };

    double B [2*2] = {
        1,2,
        1,3
    };

    double C[2*2] = {0,0,0,0};

    matrixMul(A, B, 2, 2, 2, C);

    EXPECT_DOUBLE_EQ(C[0], 3.0);
    EXPECT_DOUBLE_EQ(C[1], 8.0);
    EXPECT_DOUBLE_EQ(C[2], 7.0);
    EXPECT_DOUBLE_EQ(C[3], 18.0);

    double a[3*1] = {1,2,3};
    double b[1*2] = {1,3};
    double c[3*2];
    matrixMul(a, b, 3, 1, 2, c);
    EXPECT_DOUBLE_EQ(c[0], 1.0);
    EXPECT_DOUBLE_EQ(c[1], 3.0);
    EXPECT_DOUBLE_EQ(c[2], 2.0);
    EXPECT_DOUBLE_EQ(c[3], 6.0);
    EXPECT_DOUBLE_EQ(c[4], 3.0);
    EXPECT_DOUBLE_EQ(c[5], 9.0);
}



// To use a test fixture, derive a class from testing::Test.
class MatrixTest : public testing::Test {
 protected:  // You should make the members protected s.t. they can be
             // accessed from sub-classes.

  // virtual void SetUp() will be called before each test is run.  You
  // should define it if you need to initialize the varaibles.
  // Otherwise, this can be skipped.
  virtual void SetUp() {

  }

  // virtual void TearDown() will be called after each test is run.
  // You should define it if there is cleanup work to do.  Otherwise,
  // you don't have to provide it.
  //
  // virtual void TearDown() {
  // }

  // Declares the variables your tests want to use.
  Matrix m1;

};


TEST_F(MatrixTest, DefaultConstructor) {
	EXPECT_EQ(m1.GetRows(),0) << "The default constructor should set rows to 0!";
	EXPECT_EQ(m1.GetCols(),0) << "The default constructor should set cols to 0!";
}


TEST_F(MatrixTest, ConstructorThatCanInitialize) {
	Matrix m2(2,2);
	Matrix m4(4,4,true);
	EXPECT_EQ(m2.GetRows(),2);
	EXPECT_EQ(m2.GetCols(),2);
	EXPECT_EQ(m4.GetRows(),4);
	EXPECT_EQ(m4.GetCols(),4);
	for (int i=0; i<4; i++) {
		for (int j=0; j<4; j++) {
			EXPECT_DOUBLE_EQ(m4(i,j),0.0) << "The matrix is not initialized to 0!";
		}
	}
}

