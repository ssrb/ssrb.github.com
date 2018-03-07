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
#include "MildSlopeEquation.h"
#include "LinearSystem.h"
#include "Solver.h"
#include "Wave.h"

#include "tinyxml/tinyxml.h"

// MPI
#include <mpi.h>
#include <mpe.h>

// STL - libstdc++
#include <cstdlib>
#include <cstdio>

namespace {
	// Print device properties
	void printDevProp(cudaDeviceProp devProp)
	{
	    printf("Major revision number:         %d\n",  devProp.major);
	    printf("Minor revision number:         %d\n",  devProp.minor);
	    printf("Name:                          %s\n",  devProp.name);
	    printf("Total global memory:           %u\n",  devProp.totalGlobalMem);
	    printf("Total shared memory per block: %u\n",  devProp.sharedMemPerBlock);
	    printf("Total registers per block:     %d\n",  devProp.regsPerBlock);
	    printf("Warp size:                     %d\n",  devProp.warpSize);
	    printf("Maximum memory pitch:          %u\n",  devProp.memPitch);
	    printf("Maximum threads per block:     %d\n",  devProp.maxThreadsPerBlock);
	    for (int i = 0; i < 3; ++i)
	    printf("Maximum dimension %d of block:  %d\n", i, devProp.maxThreadsDim[i]);
	    for (int i = 0; i < 3; ++i)
	    printf("Maximum dimension %d of grid:   %d\n", i, devProp.maxGridSize[i]);
	    printf("Clock rate:                    %d\n",  devProp.clockRate);
	    printf("Total constant memory:         %u\n",  devProp.totalConstMem);
	    printf("Texture alignment:             %u\n",  devProp.textureAlignment);
	    printf("Concurrent copy and execution: %s\n",  (devProp.deviceOverlap ? "Yes" : "No"));
	    printf("Number of multiprocessors:     %d\n",  devProp.multiProcessorCount);
	    printf("Kernel execution timeout:      %s\n",  (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
	    return;
	}
	 
	int printAllDevProp()
	{
	    // Number of CUDA devices
	    int devCount;
	    cudaGetDeviceCount(&devCount);
	    printf("CUDA Device Query...\n");
	    printf("There are %d CUDA devices.\n", devCount);
	 
	    // Iterate through devices
	    for (int i = 0; i < devCount; ++i)
	    {
	        // Get device properties
	        printf("\nCUDA Device #%d\n", i);
	        cudaDeviceProp devProp;
	        cudaGetDeviceProperties(&devProp, i);
	        printDevProp(devProp);
	    }

	    size_t freeGpuMem, totalGpuMem, gpuMemlimit;
		cudaMemGetInfo(&freeGpuMem, &totalGpuMem);

		cudaDeviceGetLimit(&gpuMemlimit, cudaLimitMallocHeapSize);
		printf("Free GPU mem: %d, Total GPU mem: %d, Mem limit: %d\n", freeGpuMem, totalGpuMem, gpuMemlimit);
	}
}

int main(int ac, char **av) {

	TiXmlDocument conf( av[1] );
    bool loadOkay = conf.LoadFile();
    TiXmlElement* rootElement = conf.FirstChildElement( "babyhares" );
    
    TiXmlElement* harbourElement = rootElement->FirstChildElement( "harbour" );
    const char *geometryPath = harbourElement->FirstChildElement( "geometry" )->FirstChild()->ToText()->Value();
    const char *interfacePath = harbourElement->FirstChildElement( "interface" )->FirstChild()->ToText()->Value();
    const char *bathymetryPath = harbourElement->FirstChildElement( "bathymetry" )->FirstChild()->ToText()->Value();

	MPI_Init(&ac, &av);

	printAllDevProp();

	int mpiRank, mpiSize, mpiInfo;
	mpiInfo = MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	mpiInfo = MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

	using namespace BabyHares;

	double R = atof(harbourElement->FirstChildElement( "reflection" )->FirstChild()->ToText()->Value());
	double cf = atof(harbourElement->FirstChildElement( "friction" )->FirstChild()->ToText()->Value());

	int domainId = 1 + mpiRank;

	MeshPtr mesh(Mesh::Read(
		geometryPath, 
		interfacePath,
		bathymetryPath,
		domainId));

	Wave wave;
	TiXmlElement* waveElement = rootElement->FirstChildElement( "wave" );	
	wave.height = atof(waveElement->FirstChildElement( "height" )->FirstChild()->ToText()->Value());	
	wave.length = atof(waveElement->FirstChildElement( "length" )->FirstChild()->ToText()->Value());
	wave.direction = atof(av[2] /*waveElement->FirstChildElement( "direction" )->FirstChild()->ToText()->Value()*/) * M_PI / 180;
	
	std::vector<std::complex<double> > eta0(mesh->_vertices.size(), 0.0);
	std::complex<double> *eta0Ptr = &(eta0[0]);

	MildSlopeEquation equation(mesh);

	const int maxNonLinearLoopIter = 1;
	Solver<std::complex<double> > solver;
	for (int iter = 0; iter < maxNonLinearLoopIter; ++iter)
	{
		solver.solve(*equation.discretize(wave, R, cf, eta0Ptr), *mesh, eta0Ptr);
	}

	int rank = 2;
	std::vector<std::complex<double> > solution(mesh->_nbTotalVertices + 1); // 1 indexed


	solver.consolidateLocalSolutions(*mesh, eta0Ptr, rank);

	MPI_Finalize();

	return EXIT_SUCCESS;
}
