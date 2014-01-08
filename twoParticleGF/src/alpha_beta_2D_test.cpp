/*
 * alpha_beta_2D_test.cpp
 *
 *  Created on: Jan 7, 2014
 *      Author: pxiang
 */

#include "alpha_beta_2D.h"
#include "gtest/gtest.h"

TEST(GenerateIndexMatrix2DTest,TenByTen) {
	int xmax = 30;
	int ymax = 30;
	IntegerMatrix Index;
	QuartetListVector VtoG;
	std::vector<int> DimsOfV;

	generateIndexMatrix2D(xmax, ymax, Index, VtoG, DimsOfV);
	EXPECT_EQ(DimsOfV[5],28);
}
