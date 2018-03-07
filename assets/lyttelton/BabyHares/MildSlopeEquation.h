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
#include "ReverseConnectivityTable.h"

#include <tr1/memory>
#include <complex>

namespace BabyHares {

	struct Wave;

	template<typename K> class LinearSystem;
	typedef std::tr1::shared_ptr<LinearSystem<std::complex<double> > > MildSlopeSystemPtr;
	typedef std::tr1::shared_ptr<const LinearSystem<std::complex<double> > > MildSlopeSystemConstPtr;

	class MildSlopeEquation {
		public:
			MildSlopeEquation(const MeshConstPtr &mesh);
			MildSlopeSystemConstPtr discretize(const Wave &wave, double R, double cf, const std::complex<double> *eta0);

		private:

			typedef std::complex<double> K;

			void countNZCoefficents();
			MildSlopeSystemConstPtr computeCoefficents(const Wave &wave, double R,  double cf, const K *eta0) const;

			void walkDomainBoundary(const Wave &wave, 
									double R, 
									const double *n0s, 
									const K *ps,
									K (*v)[3][3], 
									K *load) const;
			void addOpenBoundaryEdge(	const Edge &edge, 
										const Wave &wave, 
										const double *n0s, 
										const K *ps, 
										K (*v)[3][3], 
										K *load) const;
			void addClosedBoundaryEdge(	const Edge &edge, 
										const Wave &wave, 
										double R, 
										const double *n0s, 
										const K *ps, 
										K (*v)[3][3]) const;

			int getEdgeTriangleId(const Edge &edge) const;

			ReverseConnectivityTable _rct;
			MeshConstPtr _mesh;
			int _nbNonNullCoefficentsAII, _nbNonNullCoefficentsAIG, _nbNonNullCoefficentsAGG;
	};
}
