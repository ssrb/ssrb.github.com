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

#include "Mesh.h"

#include <vector>

namespace BabyHares {
	class ReverseConnectivityTable {
		public:
			void init(const Mesh &mesh) {
				_head.resize(mesh._vertices.size(), -1);
				_next.resize(3 * mesh._triangles.size());
				for (int ti = 0, p = 0; ti < mesh._triangles.size(); ++ti) {
					for (int si = 0; si < 3; ++si, ++p) {
						int vi = mesh._triangles[ti].v[si];
						_next[p] = _head[vi];
						_head[vi] = p; // (p / 3) = ti, is the triangle number, the new head of the list of triangles for vertex vi;
					}
				}
			}

			class Iter {
				public:
					friend class ReverseConnectivityTable;
					const Iter operator++() {
						_p = _rct->_next[_p];
						return *this;
					}
					int operator*() const{
						return _p / 3;
					}
					bool operator!=(const Iter &rhs) const {
						return  _p != rhs._p || _rct != rhs._rct;
					}
				private:
					const ReverseConnectivityTable *_rct;
					int _p;
			};

			Iter find(int vi) const {	
					Iter iter;
					iter._rct = this;
					iter._p = _head[vi];
					return iter;
			}
			
			Iter end() const {
					Iter iter;
					iter._rct = this;
					iter._p = -1;
					return iter;
			}

		private:
			std::vector<int>  _head, _next;
	};
}