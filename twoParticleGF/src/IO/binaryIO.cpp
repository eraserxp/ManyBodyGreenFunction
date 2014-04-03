/*
 * binaryIO.cpp
 *
 *  Created on: Jan 2, 2014
 *      Author: pxiang
 */
#include "binaryIO.h"


// write the content of a ComplexMatrix into a binary file
// store the rows_count and column_count first
// then the content of each elements (real and imaginary parts)
int complexMatrixToBytes(ComplexMatrix& cm, std::string filename) {
	// open the file to write
	std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);
	if(!file) {
		std::cout << "Cannot open file.";
		return 1;
	}
	int rows = cm.GetRows();
	int cols = cm.GetCols();
	file.write((char*)&rows, sizeof(int));
	file.write((char*)&cols, sizeof(int));
	complex_mkl *p = cm.getPointer();
	file.write((char*)p, sizeof(complex_mkl)*rows*cols);
	file.close();
	return 0;
}


int bytesToComplexMatrix(ComplexMatrix& cm, std::string filename) {
	// open the file to write
	std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);
	if(!file) {
		std::cout << "Cannot open file.";
		return 1;
	}
	int rows;
	int cols;
	file.read((char*)&rows, sizeof(int));
	file.read((char*)&cols, sizeof(int));
	complex_mkl *p = new complex_mkl[rows*cols];
	file.read((char*)p, sizeof(complex_mkl)*rows*cols);
	file.close();
	cm = ComplexMatrix(p, rows, cols);
	delete [] p;
	return 0;
}


// save an eigen matrix into a binary file
template<typename MatrixType>
void save(std::string filename, const MatrixType& m)
{
ofstream f(filename.c_str(), ios::binary);
// write the row_count and col_count into the file
f.write((char *)&m.rows(), sizeof(m.rows()));
f.write((char *)&m.cols(), sizeof(m.cols()));
// write the matrix elements into the file
f.write((char *)&m.data(), sizeof(typename MatrixType::Scalar)*m.cols()*m.cols());
f.close();
}

// load an eigen matrix from a binary file
template<typename MatrixType>
void load(std::string filename, MatrixType& m)
{
typename MatrixType::Index rows, cols;
ifstream f(filename.c_str(), ios::binary);
f.read((char *)&rows, sizeof(rows));
f.read((char *)&cols, sizeof(cols));
m.resize(rows, cols);
f.read((char *)&m.data(), sizeof(typename MatrixType::Scalar)*rows*cols);
if (f.bad())
throw std::exception("Error reading matrix");
f.close();
}
