/*
 * binaryIO_test.cpp
 *
 *  Created on: Jan 2, 2014
 *      Author: pxiang
 */
#include "gtest/gtest.h"
#include "binaryIO.h"

TEST(ComplexMatrixIOTest,TwoByTwoMatirx) {
	ComplexMatrix cm(2,2);
	cm(0,0).real = 1.0;
	cm(0,0).imag = 0.5;
	cm(0,1).real = 1.3;
	cm(0,1).imag = 0.3;
	cm(1,0).real = 1.89;
	cm(1,0).imag = 0.35;
	cm(1,1).real = 1.07899;
	cm(1,1).imag = 0.1135;
	complexMatrixToBytes(cm,"cm.bin");
	ComplexMatrix cm2;
	bytesToComplexMatrix(cm2, "cm.bin");
	EXPECT_EQ(cm.GetRows(), cm2.GetRows());
	EXPECT_EQ(cm.GetCols(), cm2.GetCols());
	for (int i=0; i<cm.GetRows(); ++i) {
		for (int j=0; j<cm.GetCols(); ++j) {
			EXPECT_DOUBLE_EQ(cm(i,j).real, cm2(i,j).real);
			EXPECT_DOUBLE_EQ(cm(i,j).imag, cm2(i,j).imag);
		}
	}
}


TEST(ComplexMatrixIOTest,TwoByOneMatirx) {
	ComplexMatrix cm(2,1);
	cm(0,0).real = 1.0;
	cm(0,0).imag = 0.5;
	cm(1,0).real = 1.89;
	cm(1,0).imag = 0.35;
	complexMatrixToBytes(cm,"cm.bin");
	ComplexMatrix cm2;
	bytesToComplexMatrix(cm2, "cm.bin");
	EXPECT_EQ(cm.GetRows(), cm2.GetRows());
	EXPECT_EQ(cm.GetCols(), cm2.GetCols());
	for (int i=0; i<cm.GetRows(); ++i) {
		for (int j=0; j<cm.GetCols(); ++j) {
			EXPECT_DOUBLE_EQ(cm(i,j).real, cm2(i,j).real);
			EXPECT_DOUBLE_EQ(cm(i,j).imag, cm2(i,j).imag);
		}
	}
}


TEST(ComplexMatrixIOTest,OneByTwoMatirx) {
	ComplexMatrix cm(1,2);
	cm(0,0).real = 1.0;
	cm(0,0).imag = 0.5;
	cm(0,1).real = 1.3;
	cm(0,1).imag = 0.3;
	complexMatrixToBytes(cm,"cm.bin");
	ComplexMatrix cm2;
	bytesToComplexMatrix(cm2, "cm.bin");
	EXPECT_EQ(cm.GetRows(), cm2.GetRows());
	EXPECT_EQ(cm.GetCols(), cm2.GetCols());
	for (int i=0; i<cm.GetRows(); ++i) {
		for (int j=0; j<cm.GetCols(); ++j) {
			EXPECT_DOUBLE_EQ(cm(i,j).real, cm2(i,j).real);
			EXPECT_DOUBLE_EQ(cm(i,j).imag, cm2(i,j).imag);
		}
	}
}


TEST(ComplexMatrixIOTest,OneByOneMatirx) {
	ComplexMatrix cm(1,1);
	cm(0,0).real = 1.0;
	cm(0,0).imag = 0.5;
	complexMatrixToBytes(cm,"cm.bin");
	ComplexMatrix cm2;
	bytesToComplexMatrix(cm2, "cm.bin");
	EXPECT_EQ(cm.GetRows(), cm2.GetRows());
	EXPECT_EQ(cm.GetCols(), cm2.GetCols());
	for (int i=0; i<cm.GetRows(); ++i) {
		for (int j=0; j<cm.GetCols(); ++j) {
			EXPECT_DOUBLE_EQ(cm(i,j).real, cm2(i,j).real);
			EXPECT_DOUBLE_EQ(cm(i,j).imag, cm2(i,j).imag);
		}
	}
}


TEST(ComplexMatrixIOTest,ThousandByThousandMatirx) {
	ComplexMatrix cm(1000,1000);
	RandomNumberGenerator random;
	for (int i=0; i<cm.GetRows(); ++i) {
		for (int j=0; j<cm.GetCols(); ++j) {
			cm(i,j) = random.randomComplex();
		}
	}
	complexMatrixToBytes(cm,"cm.bin");
	ComplexMatrix cm2;
	bytesToComplexMatrix(cm2, "cm.bin");
	EXPECT_EQ(cm.GetRows(), cm2.GetRows());
	EXPECT_EQ(cm.GetCols(), cm2.GetCols());
	for (int i=0; i<cm.GetRows(); ++i) {
		for (int j=0; j<cm.GetCols(); ++j) {
			EXPECT_DOUBLE_EQ(cm(i,j).real, cm2(i,j).real);
			EXPECT_DOUBLE_EQ(cm(i,j).imag, cm2(i,j).imag);
		}
	}
}
