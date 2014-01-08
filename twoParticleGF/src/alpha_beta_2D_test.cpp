/*
 * alpha_beta_2D_test.cpp
 *
 *  Created on: Jan 7, 2014
 *      Author: pxiang
 */

#include "alpha_beta_2D.h"
#include "gtest/gtest.h"

TEST(GenerateIndexMatrix2DTest,SomeRandomTest) {
	int xmax = 50;
	int ymax = 50;
	IntegerMatrix Index;
	QuartetListVector VtoG;
	std::vector<int> DimsOfV;

	generateIndexMatrix2D(xmax, ymax, Index, VtoG, DimsOfV);
	// some random test
	EXPECT_EQ(DimsOfV[5],28);
	EXPECT_EQ(VtoG[5].size(),28);

	EXPECT_EQ(DimsOfV[10],140);
	EXPECT_EQ(VtoG[10].size(),140);

	EXPECT_EQ(DimsOfV[1],2);
	EXPECT_EQ(VtoG[1].size(),2);

	EXPECT_EQ(DimsOfV[50+50+50+49],2);
	EXPECT_EQ(VtoG[50+50+50+49].size(),2);

	xmax = 9;
	ymax = 9;
	generateIndexMatrix2D(xmax, ymax, Index, VtoG, DimsOfV);
	EXPECT_EQ(DimsOfV[10],138);
	EXPECT_EQ(VtoG[10].size(),138);

	QuartetList compare;
	compare.push_back(Quartet(0, 0, 0, 4));
	compare.push_back(Quartet(0, 0, 1, 3));
	compare.push_back(Quartet(0, 0, 2, 2));
	compare.push_back(Quartet(0, 0, 3, 1));
	compare.push_back(Quartet(0, 0, 4, 0));
	compare.push_back(Quartet(0, 1, 0, 3));
	compare.push_back(Quartet(0, 1, 1, 2));
	compare.push_back(Quartet(0, 1, 2, 1));
	compare.push_back(Quartet(1, 0, 0, 3));
	compare.push_back(Quartet(1, 0, 1, 2));
	compare.push_back(Quartet(1, 0, 2, 1));
	compare.push_back(Quartet(1, 0, 3, 0));
	compare.push_back(Quartet(1, 1, 0, 2));
	compare.push_back(Quartet(2, 0, 0, 2));
	compare.push_back(Quartet(2, 0, 1, 1));
	compare.push_back(Quartet(3, 0, 0, 1));

	EXPECT_EQ(VtoG[4].size(),16);
	QuartetList::iterator it;
	QuartetList::iterator it2;
	for ( it=VtoG[4].begin(), it2=compare.begin();
			it != VtoG[4].end()&&it2!=compare.end(); ++it,++it2) {
		EXPECT_TRUE(*it==*it2);
	}

}



TEST(GenerateNeighbors2DTest, ThreeByThreeCrystal) {
	/* y
	 * ^
	 * |
	 * | 6	7	8
	 * | 3	4	5
	 * | 0	1	2
	 * ------------> x
	 */
	int xmax = 2; // index starts from 0
	int ymax = 2;
	int xi, yi, xj, yj;
	Neighbors2D neighbor;
	Neighbors2D compare;


	xi = 0;
	yi = 0;
	xj = 1;
	yj = 0;
	compare = Neighbors2D(-1,-1,-1,-1,-1,2,1,-1);
	neighbor = generateNeighbors2D(xi, yi, xj, yj, xmax, ymax);
	EXPECT_TRUE(neighbor==compare);

	xi = 0;
	yi = 0;
	xj = 1;
	yj = 1;
	compare = Neighbors2D(-1,1,1,-1,0,2,2,0);
	neighbor = generateNeighbors2D(xi, yi, xj, yj, xmax, ymax);
	EXPECT_TRUE(neighbor==compare);

	xi = 0;
	yi = 1;
	xj = 1;
	yj = 1;
	compare = Neighbors2D(-1,-1,-1,0,-1,2,2,-1);
	neighbor = generateNeighbors2D(xi, yi, xj, yj, xmax, ymax);
	EXPECT_TRUE(neighbor==compare);

	xi = 1;
	yi = 1;
	xj = 2;
	yj = 2;
	compare = Neighbors2D(0,2,2,0,1,-1,-1,1);
	neighbor = generateNeighbors2D(xi, yi, xj, yj, xmax, ymax);
	EXPECT_TRUE(neighbor==compare);
}
