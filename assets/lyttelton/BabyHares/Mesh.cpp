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
#include "Mesh.h"

#include <fstream>
#include <set>

#include <tr1/unordered_set>
#include <algorithm>
#include <sstream>

#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <cmath>
#include <cassert>

 #define MeshError(FMT, ...) fprintf(stderr,"Mesh: " FMT, __VA_ARGS__);

namespace BabyHares {

	namespace {

		inline bool GotoSection(std::ifstream &stream, const char *tag) {
			std::string line;
			while(std::getline(stream, line)) {
			  	if (line == tag) {
			  		return true;
			  	}
			}
			return false;
		}

		inline bool SkipLines(std::ifstream &stream, int n) {
			int i;
			std::string skippedLine;
			for (i = 0; i < n && std::getline(stream, skippedLine); ++i);
			return i == n;
		}

		bool ReadInterfaceFile(	const char *interfaceFileName,
								const MeshPtr &mesh,
								std::vector<int> &interfaceList) {
			std::ifstream interfaceStream(interfaceFileName);
			if (!interfaceStream.is_open()) {
			 	return false;
			}
			interfaceList.reserve(1000);
			int vid;
			while (interfaceStream >> vid) {
				interfaceList.push_back(vid);
			}
			//mesh->_interface.shrink_to_fit();
			mesh->_nbInterfaceVertices = interfaceList.size();			
			return true;
		}

		bool ReadTrianglesWithinDomain(	std::ifstream &meshStream, 
										std::ifstream &depthStream,
										int domainId,
										const MeshPtr &mesh,
										const std::vector<int>  &interfaceList,
										std::set<int> &interiorVertexIds) {			
			if (!GotoSection(meshStream, "Triangles")) {
				return false;
			}
			int totalNbTriangles;
			int nbDomains = 5;
			meshStream >> totalNbTriangles;
			// Estimate the required storage based on the number of subdomains;
			std::tr1::unordered_set<int> interfaceSet(interfaceList.begin(), interfaceList.end());
			mesh->_triangles.reserve((int)round((1.1 * totalNbTriangles) / nbDomains));
			float depth;
			for (int ti = 0; ti < totalNbTriangles; ++ti) {
				Triangle triangle;
				meshStream >> triangle.v[0] >> triangle.v[1] >> triangle.v[2] >> triangle.domainId;
				depthStream >> depth;
				if (triangle.domainId == domainId) {
					mesh->_triangles.push_back(triangle);
					mesh->_depth.push_back(depth);
					for (int vi = 0; vi < 3; ++vi) {
						if (!interfaceSet.count(triangle.v[vi])) {
							interiorVertexIds.insert(triangle.v[vi]);
						}
					}
				}
			}

			printf("%u triangles\n", mesh->_triangles.size());\
			//mesh->_triangles.shrink_to_fit();
			return true;
		}

		void InitializeVertexIndexConvertionTables(	const MeshPtr &mesh,
													const std::set<int> &interiorVertexIds,
													const std::vector<int> &interfaceList) {
			mesh->_localToGlobal.reserve(interiorVertexIds.size() + interfaceList.size());
			mesh->_localToGlobal.resize(interiorVertexIds.size());
			int localId = 0;
			for (std::set<int>::const_iterator vidIter(interiorVertexIds.begin());
				vidIter != interiorVertexIds.end();
				++vidIter, ++localId) {
				mesh->_localToGlobal[localId] = *vidIter;
			}
			mesh->_localToGlobal.insert(mesh->_localToGlobal.end(), interfaceList.begin(), interfaceList.end());
			for (int i = 0; i < mesh->_localToGlobal.size(); ++i) {
				mesh->_globalToLocal.insert(std::make_pair(mesh->_localToGlobal[i], i));
			}
		}

		void RenumberTriangleVertexIndices(const MeshPtr &mesh) {
			for (int ti = 0; ti < mesh->_triangles.size(); ++ti) {
				Triangle &triangle(mesh->_triangles[ti]);
				for (int vi = 0; vi < 3; ++vi) {
					triangle.v[vi] = mesh->_globalToLocal[triangle.v[vi]];
				}
			}
		}

		bool ReadVerticesWithinDomain(std::ifstream &meshStream, const MeshPtr &mesh) {
			meshStream.seekg(0);
			if (!GotoSection(meshStream, "Vertices")) {
				return false;
			}

			std::vector<int> sortedGobalIds(mesh->_localToGlobal);
			std::sort(sortedGobalIds.begin(), sortedGobalIds.end());

			std::string dummy;

			meshStream >> mesh->_nbTotalVertices;
			std::getline(meshStream, dummy);
			assert(sortedGobalIds.size() <= mesh->_nbTotalVertices);

			mesh->_vertices.resize(sortedGobalIds.size());
			for (int vi = 0, currentLine = 1; vi < sortedGobalIds.size(); ++vi, ++currentLine) {
				int globalId = sortedGobalIds[vi];
				SkipLines(meshStream, globalId - currentLine);
				currentLine = globalId;
				Vertex &vertex(mesh->_vertices[mesh->_globalToLocal[globalId]]);
				meshStream >> vertex.x >> vertex.y >> vertex.boundaryId;
				std::getline(meshStream, dummy);
			}
			return true;
		}

		bool ReadEdgesWithinDomain(std::ifstream &meshStream, const MeshPtr &mesh) {
			meshStream.seekg(0);
			if (!GotoSection(meshStream, "Edges")) {
				return false;
			}
			std::string dummy;
			int nbEdges;
			meshStream >> nbEdges;
			std::getline(meshStream, dummy);
			mesh->_boundary.reserve(nbEdges);
			for (int ei = 0; ei < nbEdges; ++ei) {
				int v0 ,v1, boundaryId;
				meshStream >> v0 >> v1 >> boundaryId;
				if (mesh->_globalToLocal.count(v0) && mesh->_globalToLocal.count(v1)) {
					Edge edge;
					edge.v[0] = mesh->_globalToLocal[v0];
					edge.v[1] = mesh->_globalToLocal[v1];
					edge.boundaryId = boundaryId;
					mesh->_boundary.push_back(edge);
				}
			}
			return true;
		}
	}

	MeshPtr Mesh::Read(const char *meshFileName, const char *interfaceFileName, const char *depthFileName, int domainId) {

		MeshPtr mesh(new Mesh);

		std::ifstream meshStream(meshFileName);
		if (!meshStream.is_open()) {
		 	MeshError("Failed to open %s: %s\n", meshFileName, strerror(errno));
		 	return MeshPtr();
		}

		std::ifstream depthStream(depthFileName);
		if (!depthStream.is_open()) {
		 	MeshError("Failed to open %s: %s\n", depthFileName, strerror(errno));
		 	return MeshPtr();
		}

		std::vector<int>  interfaceList;
		if (!ReadInterfaceFile(interfaceFileName, mesh, interfaceList)) {
			MeshError("Failed to read interface in  %s\n", interfaceFileName);
			return  MeshPtr();
		}
		
		std::set<int> interiorVertexIds;
		if (!ReadTrianglesWithinDomain(	meshStream, 
										depthStream,
										domainId,
										mesh,
										interfaceList,
										interiorVertexIds)) {
			MeshError("Failed to read triangles in %s.\n", meshFileName);
			return  MeshPtr();
		}
		
		InitializeVertexIndexConvertionTables(mesh, interiorVertexIds, interfaceList);

		RenumberTriangleVertexIndices(mesh);

		if (!ReadVerticesWithinDomain(meshStream, mesh)) {
			MeshError("Failed to read triangles vertices in %s.\n", meshFileName);
			return MeshPtr();
		}

		if (!ReadEdgesWithinDomain(meshStream, mesh)) {
			MeshError("Failed to read edges in %s.\n", meshFileName);
			return MeshPtr();
		}

		return mesh;
	}
}
