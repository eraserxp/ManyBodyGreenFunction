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
	Eigen::SelfAdjointEigenSolver<Matrix2f> eigensolver(hamiltonian);
	if (eigensolver.info() == Eigen::Success) {
		eigenValues = eigensolver.eigenvalues();
		eigenVectors = eigensolver.eigenvectors();
	} else {
		printf("Failed to solve the eigen problem!");
	}
}
