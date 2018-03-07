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
#include "Solver.h"

namespace BabyHares {

	template<>
	void Solver<std::complex<double> >::conjugateGradient(const SolverContext &ctx, const MatrixProductFunctor_t &A, const std::vector<K> &b, K *x, double tol, int maxIter) const {
		const size_t n = b.size();
		std::fill_n(x, n, 0);
		std::vector<K> r(b);
		std::vector<K> p(r);
		K rsold;
		cblas_zdotu_sub(n, &r[0], 1, &r[0], 1, &rsold);
		std::vector<K> Ap(n, 0.0);
		for (int iter = 0; iter < maxIter; ++iter) {
			A(ctx, p, &Ap);
			
			K pAp;
			cblas_zdotu_sub(n, &Ap[0], 1, &p[0], 1, &pAp);
			
			K alpha = rsold / pAp;
			cblas_zaxpy(n, &alpha, &p[0], 1, x, 1);
			
			K minusAlpha = -alpha;
			cblas_zaxpy(n, &minusAlpha, &Ap[0], 1, &r[0], 1);
			
			K rsnew;
			cblas_zdotu_sub(n, &r[0], 1, &r[0], 1, &rsnew);
			if (sqrt(std::abs(rsnew)) < tol) {
				printf("Converged after %d iterations\n", iter);
				break;
			}

			K scal = rsnew / rsold;
			cblas_zscal (n, &scal, &p[0], 1);

			K one = 1.0;
			cblas_zaxpy(n, &one, &r[0], 1, &p[0], 1);
			
			rsold = rsnew;
		}
	}

}
