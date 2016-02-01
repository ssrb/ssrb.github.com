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

#include "MeshPtr.h"

#include <cmath>
#include <vector>
#include <tr1/unordered_map>

namespace BabyHares {

	enum BoundaryId {
		kClosedBoundaryId = 1,
		kOpenBoundaryId = 2
	};

	struct Triangle {
		int v[3], domainId;
	};

	struct Vertex {
		union {
			double coord[2];
			struct {
				double x, y;
			};
		};
		int boundaryId;
	};

	struct Edge {
		int v[2], boundaryId;
	};
		
	

	class Mesh {
		public:
			static MeshPtr Read(const char *meshFileName, const char *interfaceFileName, const char *depthFileName, int domainId);
			
			size_t NbLocalInterior() const
			{
				return _vertices.size() - _nbInterfaceVertices;
			}

			size_t NbGlobalInterior() const
			{
				return _nbTotalVertices - _nbInterfaceVertices;
			}

		//private:


			std::vector<Triangle> _triangles;
			std::vector<Vertex> _vertices;
			std::vector<Edge> _boundary;
			std::vector<float> _depth;

			std::vector<int> _localToGlobal;
			std::tr1::unordered_map<int, int> _globalToLocal;
			int _nbInterfaceVertices;
			int _nbTotalVertices;

	};
}