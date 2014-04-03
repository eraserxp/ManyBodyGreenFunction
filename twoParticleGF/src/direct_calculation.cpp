/*
 * direct_calculation.cpp
 *
 *  Created on: Jan 10, 2014
 *      Author: pxiang
 */

// calculate the Green's function directly from eigenfunctions

// forming all the basis sets
#include "direct_calculation.h"


// nmax is the max index for a lattice site (index starting from 0)
// IndexMatrix(i,j) = nth means that the basis set (i,j) is the nth basis set
void formBasisSets(int nmax, IMatrix& indexMatrix, PairVector& basisSets) {
	// the total number of site
	int nsite = nmax + 1;
	// figure out the number of basis sets
	int num_basis = nmax*(1 + nmax)/2;
	indexMatrix = IMatrix(nsite, nsite);
	// allocate the space needed to store all basis sets
	basisSets.reserve(num_basis);

	// generate all basis sets
	// to avoid double counting, we only consider (n1, n2) where n1<n2
	int nth = 0;
	for (int n1=0; n1<=nmax-1; ++n1) {
		for (int n2=n1+1; n2<=nmax; ++n2) {
			indexMatrix(n1, n2) = nth;
			nth++;
			basisSets.push_back(Pair(n1, n2));
		}
	}
}


void formHamiltonianMatrix(InputParameters& pars, DMatrix& hamiltonian, IMatrix& indexMatrix,
		                   PairVector& basisSets) {
	int size = basisSets.size();
	hamiltonian = DMatrix::Zero(size, size);
	// fill in the matrix column by column  (by default, the storage order is column-major)
	for (int col=0; col<size; ++col) {
		Pair ket = basisSets[col];
		int ket_n1 = ket.First();
		int ket_n2 = ket.Second();
		// fill in the diagonal part
		if (ket_n2-ket_n1==1) {
			hamiltonian(col, col)+=pars.d0;
		}
		// fill in the off-diagonal part
		int bra_n1;
		int bra_n2;
		int row;
		// find out which basis set will produce non-interacting matrix with ket
		Neighbors ns = generateNeighbors(ket_n1, ket_n2, pars.nmax);
		if (ns.firstLeft!=-1) { // the first particle has left neighbor
			bra_n1 = ket_n1 - 1;
			bra_n2 = ket_n2;
			row = indexMatrix(bra_n1, bra_n2);
			hamiltonian(row, col) = pars.t0;
		}

		if (ns.firstRight!=-1) { // the first particle has right neighbor
			bra_n1 = ket_n1 + 1;
			bra_n2 = ket_n2;
			row = indexMatrix(bra_n1, bra_n2);
			hamiltonian(row, col) = pars.t0;
		}

		if (ns.secondLeft!=-1) { // the second particle has left neighbor
			bra_n1 = ket_n1;
			bra_n2 = ket_n2 - 1;
			row = indexMatrix(bra_n1, bra_n2);
			hamiltonian(row, col) = pars.t0;
		}

		if (ns.secondRight!=-1) { // the second particle has right neighbor
			bra_n1 = ket_n1;
			bra_n2 = ket_n2 + 1;
			row = indexMatrix(bra_n1, bra_n2);
			hamiltonian(row, col) = pars.t0;
		}

	}

	// adding the impurities
	for (int i=0; i<pars.impurityIndex.size(); ++i) {
		int nth = indexMatrix(i,i);
		hamiltonian(nth, nth) = pars.impurityStrength[i];
	}
}


void obtainEigenVectors(DMatrix& hamiltonian, DVector& eigenValues, DMatrix& eigenVectors) {
	Eigen::SelfAdjointEigenSolver<DMatrix> eigensolver(hamiltonian);
	if (eigensolver.info() == Eigen::Success) {
		eigenValues = eigensolver.eigenvalues();
		eigenVectors = eigensolver.eigenvectors();
	} else {
		printf("Failed to solve the eigen problem!");
	}
}


//helper functions for calculating the green function G(n1f, n2f, n1i, n2i, z)
void numeratorHelper(int n1f, int n2f, int n1i, int n2i, IMatrix& indexMatrix,
		             DMatrix& eigenVectors, CDArray& numerator) {
	int bra_index = indexMatrix(n1f, n2f);
	int ket_index = indexMatrix(n1i, n2i);
	DVector wavefunc_f = eigenVectors.row(bra_index);
	DVector wavefunc_i = eigenVectors.row(ket_index);
	DArray tmp = wavefunc_f.array()*wavefunc_i.array();
	numerator = tmp.cast< dcomplex >();
}

void denominatorHelper(dcomplex z, DVector& eigenValues, CDArray& oneOverDenominator) {
	int n = eigenValues.rows();
	oneOverDenominator = CDArray(n);
	for (int i=0; i<n; ++i) {
		oneOverDenominator(i) = 1.0/(z - eigenValues(i));
	}

}

dcomplex greenFunc(CDArray& numerator, CDArray& oneOverDenominator) {
	CDArray tmp = numerator*oneOverDenominator;
	return tmp.sum();
}



// calculate the density of state at (ni1, ni2)
void densityOfState_direct(InputParameters& pars, int n1i, int n2i, std::vector<dcomplex >& zList,
		                   std::vector<double>& dosList) {
	int nmax = pars.nmax;
	IMatrix indexMatrix;
	PairVector basisSets;
	formBasisSets(nmax, indexMatrix, basisSets);

	DMatrix hamiltonian;
	formHamiltonianMatrix(pars, hamiltonian, indexMatrix, basisSets);

	DVector eigenValues;
	DMatrix eigenVectors;
	obtainEigenVectors(hamiltonian, eigenValues, eigenVectors);

	int n1f = n1i;
	int n2f = n2i;
	CDArray numerator;
	numeratorHelper(n1f, n2f, n1i, n2i, indexMatrix, eigenVectors, numerator);

	dosList.clear();
	dosList.reserve(zList.size());
	for (int i=0; i<zList.size(); ++i) {
		dcomplex z = zList[i];
		CDArray oneOverDenominator;
		denominatorHelper(z, eigenValues, oneOverDenominator);
		dcomplex gf = greenFunc(numerator, oneOverDenominator);
		dosList[i] = -gf.imag()/M_PI;
	}

}
