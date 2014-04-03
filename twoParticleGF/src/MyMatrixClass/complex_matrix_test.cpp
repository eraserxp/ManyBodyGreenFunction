/*
 * complex_matrix_test.cpp
 *
 *  Created on: Dec 25, 2013
 *      Author: pxiang
 */

#include "matrix.h"
#include "complex_matrix.h"
#include "gtest/gtest.h"

TEST(MinusOperatorTest, RealMinusComplex) {
	Matrix I = Identity(2);
	ComplexMatrix cm(2,2);
	cm(0,0).real = 1.0;
	cm(0,0).imag = 0.1;
	cm(0,1).real = 2.0;
	cm(0,1).imag = 0.2;
	cm(1,0).real = 3.0;
	cm(1,0).imag = 0.3;
	cm(1,1).real = 4.0;
	cm(1,1).imag = 0.4;

	ComplexMatrix result = I - cm;
	EXPECT_DOUBLE_EQ(result(0,0).real, 0.0);
	EXPECT_DOUBLE_EQ(result(0,0).imag, -0.1);
	EXPECT_DOUBLE_EQ(result(0,1).real, -2.0);
	EXPECT_DOUBLE_EQ(result(0,1).imag, -0.2);
	EXPECT_DOUBLE_EQ(result(1,0).real, -3.0);
	EXPECT_DOUBLE_EQ(result(1,0).imag, -0.3);
	EXPECT_DOUBLE_EQ(result(1,1).real, -3.0);
	EXPECT_DOUBLE_EQ(result(1,1).imag, -0.4);
}
