/*
 * binaryIO.h
 *
 *  Created on: Jan 2, 2014
 *      Author: pxiang
 */

#ifndef BINARYIO_H_
#define BINARYIO_H_

#include <iostream>
#include <fstream>
#include <string>
#include "complex_matrix.h"
#include "random_generator.h"

int complexMatrixToBytes(ComplexMatrix& cm, std::string filename);
int bytesToComplexMatrix(ComplexMatrix& cm, std::string filename);

#endif /* BINARYIO_H_ */
