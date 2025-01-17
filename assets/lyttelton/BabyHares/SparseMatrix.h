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

#include <vector>

namespace BabyHares {
	// CSR format, index base 0
	template <typename K>
	struct SparseMatrix {

		enum Type {
			TypeGeneral,
			TypeSymmetric,
			TypeHermitian
		};

		SparseMatrix() {
			reset(TypeGeneral, 0, 0, 0);
		}

		SparseMatrix(Type type, int nrow, int ncol, int nnz) {
			reset(type, nrow, ncol, nnz);
		}

		void reset(Type iType, int iNrow, int iNcol, int iNnz) {
			type = iType;
			nrow = iNrow;
			ncol = iNcol;
			rowPtr.resize(nrow + 1);
			rowPtr[0] = 0;
			col.resize(iNnz, 0);
			val.resize(iNnz, 0.0);
		}

		Type type;
		int nrow, ncol;
		std::vector<int> rowPtr, col;
		std::vector<K> val;
	};
}