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

// save an eigen matrix into a binary file
template<typename MatrixType>
void save(std::string filename, const MatrixType& m);

// load an eigen matrix from a binary file
template<typename MatrixType>
void load(std::string filename, MatrixType& m);

#endif /* BINARYIO_H_ */
