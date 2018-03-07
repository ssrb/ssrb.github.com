// Copyright (c) 2015, Sebastien Sydney Robert Bigot
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The views and conclusions contained in the software and documentation are those
// of the authors and should not be interpreted as representing official policies,
// either expressed or implied, of the FreeBSD Project.
#include "LU.h"

namespace BabyHares {
	LUConstPtr LU::Factorize(const SparseMatrix<std::complex<double> > &mat) {

		LUPtr lu(new LU);

		lu->_mtype = 6; // complex and symmetric

		pardisoinit(lu->_pardisoHandle, &lu->_mtype, lu->_iparm);

		//pardisoLU.iparm[26] = 1; // Check input
		lu->_iparm[27] = 0; // Pardiso double precision
		lu->_iparm[34] = 1; // 0 based indexing (C-style)

		// Intel doc : Maximum number of factors with identical nonzero sparsity structure 
		// that must be keep at the same time in memory. 
		// In most applications this value is equal to 1.
		MKL_INT maxfct = 1;

		// Intel doc: Indicates the actual matrix for the solution phase. 
		// With this scalar you can define which matrix to factorize.
		// In most applications this value is 1.
		MKL_INT mnum = 1;

		MKL_INT phase = 12; // Analysis, numerical factorization

		lu->_n = mat.nrow;
		
		// Intel doc: Holds the permutation vector of size n. 
		// You can use it to apply your own fill-in reducing ordering to the solver. 
		// The permutation vector perm is used by the solver if iparm(5) = 1.
		MKL_INT *perm = 0;

		// Number of right-hand sides that need to be solved for.
		MKL_INT nrhs = 0;

		MKL_INT msglvl = 0; // Debug output => 1

		MKL_INT error = 0;

		lu->_val = &mat.val[0];
		lu->_rowPtr = &mat.rowPtr[0];
		lu->_col = &mat.col[0];

		pardiso(lu->_pardisoHandle,
				&maxfct,
				&mnum,
				&lu->_mtype,
				&phase,
				&lu->_n,
				const_cast<double *>(reinterpret_cast<const double *>(&mat.val[0])),
				const_cast<int *>(&mat.rowPtr[0]),
				const_cast<int *>(&mat.col[0]),
				perm,
				&nrhs,
				lu->_iparm,
				&msglvl,
				0,
				0,
				&error);

		return lu;
	}

	LU::~LU() {
		// Intel doc : Maximum number of factors with identical nonzero sparsity structure 
		// that must be keep at the same time in memory. 
		// In most applications this value is equal to 1.
		MKL_INT maxfct = 1;

		// Intel doc: Indicates the actual matrix for the solution phase. 
		// With this scalar you can define which matrix to factorize.
		// In most applications this value is 1.
		MKL_INT mnum = 1;

		MKL_INT phase = -1; // Release all internal memory for all matrices

		// Intel doc: Holds the permutation vector of size n. 
		// You can use it to apply your own fill-in reducing ordering to the solver. 
		// The permutation vector perm is used by the solver if iparm(5) = 1.
		MKL_INT *perm = 0;

		// Number of right-hand sides that need to be solved for.
		MKL_INT nrhs = 0;

		MKL_INT msglvl = 0; // Debug output => 1

		MKL_INT error = 0;
		pardiso(const_cast<void **>(_pardisoHandle),
				&maxfct,
				&mnum,
				const_cast<MKL_INT *>(&_mtype),
				&phase, 
				const_cast<MKL_INT *>(&_n),
				const_cast<double *>(reinterpret_cast<const double *>(_val)),
				const_cast<int *>(_rowPtr),
				const_cast<int *>(_col),
				perm,
				&nrhs,
				const_cast<MKL_INT *>(_iparm),
				&msglvl,
				0,
				0,
				&error);
	}

	void LU::solve(std::complex<double> *bx) const
	{
		// Intel doc : Maximum number of factors with identical nonzero sparsity structure 
		// that must be keep at the same time in memory. 
		// In most applications this value is equal to 1.
		MKL_INT maxfct = 1;

		// Intel doc: Indicates the actual matrix for the solution phase. 
		// With this scalar you can define which matrix to factorize.
		// In most applications this value is 1.
		MKL_INT mnum = 1;

		MKL_INT phase = 33; // Solve

		// Intel doc: Holds the permutation vector of size n. 
		// You can use it to apply your own fill-in reducing ordering to the solver. 
		// The permutation vector perm is used by the solver if iparm(5) = 1.
		MKL_INT *perm = 0;

		// Number of right-hand sides that need to be solved for.
		MKL_INT nrhs = 1;

		MKL_INT msglvl = 0; // Debug output

		MKL_INT error = 0;

		_iparm[5] = 1; // in-place

		std::vector<std::complex<double> > x(_n);

		pardiso(const_cast<void **>(_pardisoHandle),
				&maxfct,
				&mnum,
				const_cast<MKL_INT *>(&_mtype),
				&phase, 
				const_cast<MKL_INT *>(&_n),
				const_cast<double *>(reinterpret_cast<const double *>(_val)),
				const_cast<int *>(_rowPtr),
				const_cast<int *>(_col),
				perm,
				&nrhs,
				_iparm,
				&msglvl,
				(void *)bx,
				(void *)&x[0],
				&error);
	}

	void LU::solve(const std::vector<std::complex<double> > &b, std::complex<double> *x) const {
		// Intel doc : Maximum number of factors with identical nonzero sparsity structure 
		// that must be keep at the same time in memory. 
		// In most applications this value is equal to 1.
		MKL_INT maxfct = 1;

		// Intel doc: Indicates the actual matrix for the solution phase. 
		// With this scalar you can define which matrix to factorize.
		// In most applications this value is 1.
		MKL_INT mnum = 1;

		MKL_INT phase = 33; // Solve

		// Intel doc: Holds the permutation vector of size n. 
		// You can use it to apply your own fill-in reducing ordering to the solver. 
		// The permutation vector perm is used by the solver if iparm(5) = 1.
		MKL_INT *perm = 0;

		// Number of right-hand sides that need to be solved for.
		MKL_INT nrhs = 1;

		MKL_INT msglvl = 0; // Debug output

		MKL_INT error = 0;

		_iparm[5] = 0;

		pardiso(const_cast<void **>(_pardisoHandle),
				&maxfct,
				&mnum,
				const_cast<MKL_INT *>(&_mtype),
				&phase, 
				const_cast<MKL_INT *>(&_n),
				const_cast<double *>(reinterpret_cast<const double *>(_val)),
				const_cast<int *>(_rowPtr),
				const_cast<int *>(_col),
				perm,
				&nrhs,
				_iparm,
				&msglvl,
				(void *)&b[0],
				x,
				&error);
	}
}
