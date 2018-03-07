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
#pragma once

// BabyHares
#include "SparseMatrix.h"

#include <complex>

// cuda
#include "cusparse_v2.h"
#include <thrust/device_vector.h>

namespace BabyHares {
	// CSR format, index base 0
	template<typename K>
	struct CudaSparseMatrix {

		typedef K K;

		CudaSparseMatrix(const SparseMatrix<K> &matHost) :
			nrow(matHost.nrow),
			ncol(matHost.ncol),
			nval(matHost.val.size()),
			rowPtr(matHost.rowPtr.begin(), matHost.rowPtr.end()),
			col(matHost.col.begin(), matHost.col.end()) {

			cudaMalloc((void **)&val, nval * sizeof(K));
			cudaMemcpy(val, &matHost.val[0], nval * sizeof(K), cudaMemcpyHostToDevice);
			
		    cusparseCreateMatDescr(&descr);
		    if (matHost.type == SparseMatrix<K>::TypeSymmetric) {
		    	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_SYMMETRIC);
		    	cusparseSetMatFillMode(descr, CUSPARSE_FILL_MODE_UPPER);
	            cusparseSetMatDiagType(descr, CUSPARSE_DIAG_TYPE_NON_UNIT);
		    } else {
		    	cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
		    }
		    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
		}

		virtual ~CudaSparseMatrix() {
		    cusparseDestroyMatDescr(descr);
		    cudaFree(val);
		}

		void mv(cusparseHandle_t cusparseHandle,
				cusparseOperation_t trans,
				const K &alpha,
				const K *xdevice,
				const K &beta,
				K *ydevice) const;

		int nrow, ncol, nval;
		cusparseMatDescr_t descr;
	    thrust::device_vector<int> rowPtr;
		thrust::device_vector<int> col;
		K *val;
	};

	template<>
	void CudaSparseMatrix<std::complex<double> >::mv(cusparseHandle_t cusparseHandle, cusparseOperation_t trans, const K &alpha, const K *xdevice, const K &beta, K *ydevice) const {
			cusparseStatus_t status = cusparseZcsrmv(	cusparseHandle, 
														trans, 
														nrow,
														ncol,
														nval,
														(cuDoubleComplex *)&alpha, 
														descr,
														(cuDoubleComplex *)val,
														thrust::raw_pointer_cast(rowPtr.data()), 
														thrust::raw_pointer_cast(col.data()),
														(cuDoubleComplex *)xdevice, 
														(cuDoubleComplex *)&beta, 
														(cuDoubleComplex *)ydevice);
	}


}