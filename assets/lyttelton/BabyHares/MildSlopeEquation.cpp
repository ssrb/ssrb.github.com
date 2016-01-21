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
#include "MildSlopeEquation.h"
#include "Mesh.h"
#include "Wave.h"
#include "SparseMatrix.h"
#include "LinearSystem.h"
#include "GpuAssembly.h"

#include "mkl.h"

#include <cassert>

#include <vector>
#include <algorithm>
#include <tr1/unordered_set>

namespace BabyHares {

typedef ReverseConnectivityTable::Iter RCTI;

namespace {
	template<typename T>
	inline void CommitSparseMatrixRow(SparseMatrix<T> &mat, int row, int nnz, int colOffset, const std::vector<T> &rowVal) {
		mat.rowPtr[row + 1] = nnz;
		std::sort(&mat.col[mat.rowPtr[row]], &mat.col[mat.rowPtr[row + 1]]);
		for (int k = mat.rowPtr[row]; k < mat.rowPtr[row + 1]; ++k) {
			mat.val[k] = rowVal[mat.col[k]];
			mat.col[k] -= colOffset;
		}
	}

	inline int GetTriangleVertexId(const Triangle &triangle, int meshVertexId) {
		return meshVertexId == triangle.v[0] ? 0 :
						(meshVertexId == triangle.v[1] ? 1 : 2);
	}
}

MildSlopeEquation::MildSlopeEquation(const MeshConstPtr &mesh)
: _mesh(mesh)
{
	// Initialize a reverse connectivity table 
	// used to process the linear system one row at a time
	_rct.init(*_mesh);

	// Perform a first pass on the mesh to count 
	// the non-zero coefficient of the linear system
	countNZCoefficents();
}

MildSlopeSystemConstPtr MildSlopeEquation::discretize(const Wave &wave, double R, double cf, const std::complex<double> *eta0) {
	
	// We know how much memory neeeds to be allocated now
	return computeCoefficents(wave, R, cf, eta0);
}

void  MildSlopeEquation::countNZCoefficents() {
	
	int nbInteriorVertices = _mesh->_vertices.size() - _mesh->_nbInterfaceVertices;

	_nbNonNullCoefficentsAII = 0; 
	_nbNonNullCoefficentsAIG = 0;
	_nbNonNullCoefficentsAGG = 0;

	std::vector<int> color(_mesh->_vertices.size(), -1);
	for (int row = 0; row < _mesh->_vertices.size(); ++row) {
		const Vertex &vi(_mesh->_vertices[row]);
		// Inspect each triangle connected to vertex row
		for (RCTI ti = _rct.find(row); ti != _rct.end(); ++ti) {
			const Triangle &triangle(_mesh->_triangles[*ti]);
			for (int sj = 0; sj < 3; ++sj) {
				int column = triangle.v[sj];
				const Vertex &vj(_mesh->_vertices[column]);
				// Do not count coefficient twice
				if (color[column] != row) {
					color[column] = row;
					// AII or AIG
					if (row < nbInteriorVertices) {			
							// AII			
							if (column < nbInteriorVertices) {
								// AII is symmetric
								if (row <= column) {
									++_nbNonNullCoefficentsAII;
								}
							// AIG
							} else {
								++_nbNonNullCoefficentsAIG;
							}
					// AGG
					} else {
						// AGG is symmetric
						if (row <= column) {
							++_nbNonNullCoefficentsAGG;
						}
					}
				}
			}
		}
	}
}

void MildSlopeEquation::walkDomainBoundary(
	const Wave &wave, 
	double R, 
	const double *n0s, 
	const K *ps,
	K (*v)[3][3], 
	K *load
) const {
	for (int ei = 0; ei < _mesh->_boundary.size(); ++ei) {
		const Edge &edge(_mesh->_boundary[ei]);
		switch (edge.boundaryId) {
			case kOpenBoundaryId:
				addOpenBoundaryEdge(edge, wave, n0s, ps, v, load);
				break;
			case kClosedBoundaryId:
				addClosedBoundaryEdge(edge, wave, R, n0s, ps, v);
				break;
		}
	}
}

void MildSlopeEquation::addOpenBoundaryEdge(
	const Edge &edge, 
	const Wave &wave, 
	const double *n0s, 
	const K *ps, 
	K (*v)[3][3], 
	K *load
) const {

	int triangleIndex = getEdgeTriangleId(edge),
		si = GetTriangleVertexId(_mesh->_triangles[triangleIndex], edge.v[0]),
		sj = GetTriangleVertexId(_mesh->_triangles[triangleIndex], edge.v[1]);

	double k0 =  wave.waveNumber();

	const Vertex &vi(_mesh->_vertices[edge.v[0]]);
	const Vertex &vj(_mesh->_vertices[edge.v[1]]);
	double dx = vj.coord[0] - vi.coord[0], dy = vj.coord[1] - vi.coord[1];
	double len = sqrt(dx * dx + dy * dy);
	double n0 =  n0s[triangleIndex];	
	K p = ps[triangleIndex];

	K prefactor(0, -1.0 / (k0 * k0));

	v[triangleIndex][si][sj] += prefactor * n0 * (p * len / 6 + 1 / (2 * p * len));
	v[triangleIndex][sj][si] += prefactor * n0 * (p * len / 6 + 1 / (2 * p * len));
	v[triangleIndex][si][si] += prefactor * n0 * (p * len / 3 - 1 / (2 * p * len));
	v[triangleIndex][sj][sj] += prefactor * n0 * (p * len / 3 - 1 / (2 * p * len));

	K wavei = wave.elevationAt(vi.x, vi.y);
	K wavej = wave.elevationAt(vj.x, vj.y);
	
	load[edge.v[0]] += -prefactor * n0 * (wavej - wavei) / (2 * p * len);
	load[edge.v[1]] += -prefactor * n0 * (wavei - wavej) / (2 * p * len);

	float xwavedirection = (dy * cos(wave.direction) - dx * sin(wave.direction)) / len;

	load[edge.v[0]] += -prefactor * n0 * p * (xwavedirection - 1) * (len / 6) * (2 * wavei + wavej);
	load[edge.v[1]] += -prefactor * n0 * p * (xwavedirection - 1) * (len / 6) * (wavei + 2 * wavej);
	
}

void MildSlopeEquation::addClosedBoundaryEdge(
	const Edge &edge, 
	const Wave &wave, 
	double R, 
	const double *n0s, 
	const K *ps, 
	K (*v)[3][3]
) const {

	int triangleIndex = getEdgeTriangleId(edge),
		si = GetTriangleVertexId(_mesh->_triangles[triangleIndex], edge.v[0]),
		sj = GetTriangleVertexId(_mesh->_triangles[triangleIndex], edge.v[1]);

	double k0 =  wave.waveNumber();

	const Vertex &vi(_mesh->_vertices[edge.v[0]]);
	const Vertex &vj(_mesh->_vertices[edge.v[1]]);
	double dx = vj.coord[0] - vi.coord[0], dy = vj.coord[1] - vi.coord[1];
	double len = sqrt(dx * dx + dy * dy);
	double n0 =  n0s[triangleIndex];
	K p = ps[triangleIndex];

	K prefactor(0, -(1.0 - R) / (k0 * k0 * (1.0 + R)));

	v[triangleIndex][si][sj] += prefactor * n0 * (p * len / 6 + 1 / (2 * p * len));
	v[triangleIndex][sj][si] += prefactor * n0 * (p * len / 6 + 1 / (2 * p * len));
	v[triangleIndex][si][si] += prefactor * n0 * (p * len / 3 - 1 / (2 * p * len));
	v[triangleIndex][sj][sj] += prefactor * n0 * (p * len / 3 - 1 / (2 * p * len));
}

int MildSlopeEquation::getEdgeTriangleId(const Edge &edge) const {
	std::tr1::unordered_set<int> tids;
	for (RCTI ti = _rct.find(edge.v[0]); ti != _rct.end(); ++ti) {
		tids.insert(*ti);
	}
	for (RCTI ti = _rct.find(edge.v[1]); ti != _rct.end(); ++ti) {
		if (tids.count(*ti)) {
			return *ti;
		}
	}
	assert(false);
}

MildSlopeSystemConstPtr  MildSlopeEquation::computeCoefficents(const Wave &wave, double R, double cf, const K *eta0) const {

	MildSlopeSystemPtr ls(new LinearSystem<K>);

	int coeffAII = 0, coeffAIG = 0, coeffAGG = 0;
	std::vector<int> color(_mesh->_vertices.size(), -1);
	std::vector<K> coefficients(_mesh->_vertices.size(), 0.0);

	int nbInteriorVertices = _mesh->_vertices.size() - _mesh->_nbInterfaceVertices;

	ls->AII.reset(SparseMatrix<K>::TypeSymmetric, nbInteriorVertices, nbInteriorVertices, _nbNonNullCoefficentsAII);
	ls->AIG.reset(SparseMatrix<K>::TypeGeneral, nbInteriorVertices, _mesh->_nbInterfaceVertices, _nbNonNullCoefficentsAIG);
	ls->AGG.reset(SparseMatrix<K>::TypeSymmetric, _mesh->_nbInterfaceVertices, _mesh->_nbInterfaceVertices, _nbNonNullCoefficentsAGG);
	ls->loadVector.resize(_mesh->_vertices.size(), 0.0);

	std::vector<K> v(9 * _mesh->_triangles.size());
	std::vector<double> n0s(_mesh->_triangles.size());
	std::vector<K> ps(_mesh->_triangles.size());
	
	ComputeInteriorCoefficientsOnGpu(wave, *_mesh, cf, eta0, v, n0s, ps);
	
	K (*vv)[3][3] = (K (*)[3][3])(&v[0]);
	walkDomainBoundary(wave, R, &(n0s[0]), &(ps[0]), vv, &(ls->loadVector[0]));

	for (int row = 0; row < _mesh->_vertices.size(); ++row) {
		const Vertex &vi(_mesh->_vertices[row]);
		for (RCTI ti = _rct.find(row); ti != _rct.end(); ++ti) {

			int triangleIndex = *ti;

			const Triangle &triangle(_mesh->_triangles[triangleIndex]);

			int si = row == triangle.v[0] ? 0 :
					(row == triangle.v[1] ? 1 : 2);

			for (int sj = 0; sj < 3; ++sj) {
				int column = triangle.v[sj];
				const Vertex &vj(_mesh->_vertices[column]);
				
				// But insert once
				if (color[column] != row) {
					color[column] = row;
					coefficients[column] = 0.0;
					// AII or AIG
					if (row < nbInteriorVertices) {
						// AII
						if (column < nbInteriorVertices) {
							if (row <= column) {
								ls->AII.col[coeffAII++] = column;
							}
						// AIG
						} else {
							ls->AIG.col[coeffAIG++] = column;
						}						
					} else  {
						// AGG is symmetrical
						if (row <= column) {
							ls->AGG.col[coeffAGG++] = column;
						}
					}
				}

				// Accumulate contributions
				coefficients[column] += vv[triangleIndex][si][sj];
			}
		}

		if (row < nbInteriorVertices) {
			CommitSparseMatrixRow(ls->AII, row, coeffAII, 0, coefficients);
			CommitSparseMatrixRow(ls->AIG, row, coeffAIG, nbInteriorVertices, coefficients);
		} else {
			CommitSparseMatrixRow(ls->AGG, row - nbInteriorVertices, coeffAGG, nbInteriorVertices, coefficients);
		}
	}

	return ls;
}

}
