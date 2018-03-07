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
#include "LinearSystem.h"
#include "Mesh.h"
#include "LU.h"
#include "CudaSparseMatrix.h"

// MPI
#include <mpi.h>

// Pardiso + vectorized math
#include "mkl.h"

// cuSPARSE
#include <cuda_runtime.h>
#include <thrust/host_vector.h>

// STL -libstdc++
#include <numeric>
#include <tr1/functional>
#include <tr1/memory>
#include <fstream>

namespace
{
	void DumpVector(const double *x, size_t n, const std::string &name)
	{
		int mpiRank;
		MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
		std::stringstream myName;
		myName << name << "." << mpiRank;
		std::ofstream out(myName.str().c_str());
		for (int i = 0; i < (int)n; ++i)
		{
			out << x[i] << std::endl;
		}
	}
}

namespace BabyHares {

	struct SolverContext;

	template<typename K>
	class Solver {
		public:
			// Solve _localy_
			void solve(const LinearSystem<K> &ls, const Mesh &mesh, K *x) const;
			void consolidateLocalSolutions(const Mesh &mesh, const K *localSol, int rank) const;
			typedef K K;
		private:

			struct SolverContext {

				SolverContext(const LinearSystem<K> &linearSystem)	: ls(linearSystem) {
					int workBuffSize = std::max(ls.AIG.nrow, ls.AIG.ncol) * sizeof(K);
					cusparseCreate(&cusparseHandle);
					K *tmp;
					cudaMalloc(&tmp, workBuffSize);
					workGpuBuff1 = thrust::device_ptr<K>(tmp);
					cudaMalloc(&tmp, workBuffSize);
					workGpuBuff2 = thrust::device_ptr<K>(tmp);
					cudaMalloc(&tmp, workBuffSize);
					workGpuBuff3 = thrust::device_ptr<K>(tmp);
				}

				virtual ~SolverContext() {
					cusparseDestroy(cusparseHandle);
					cudaFree(workGpuBuff1.get());
					cudaFree(workGpuBuff2.get());
					cudaFree(workGpuBuff3.get());
				}

				const LinearSystem<K> &ls;
				LUConstPtr factorizedIStiffness;

				typedef std::tr1::shared_ptr<const CudaSparseMatrix<K> > CudaSparseMatrixConstPtr;

				CudaSparseMatrixConstPtr AIGDevice;
				CudaSparseMatrixConstPtr AGGDevice;
				cusparseHandle_t cusparseHandle;

				thrust::device_ptr<K> workGpuBuff1, workGpuBuff2, workGpuBuff3;
			};

			void computeBTilde(const SolverContext &ctx, std::vector<K> *bt) const;

			typedef std::tr1::function<void (const SolverContext &ctx, const std::vector<K> &p, std::vector<K> *Ap)> MatrixProductFunctor_t;
			void conjugateGradient(const SolverContext &ctx, const MatrixProductFunctor_t &A, const std::vector<K> &b, K *x, double tol, int maxIter) const;
			void multiplyBySchurComplement(const SolverContext &ctx, const std::vector<K> &p, std::vector<K> *Ap) const;

			void solveLocalInterior(const SolverContext &ctx, K *x) const;
	};
}

namespace BabyHares {

	template<typename K>
	void Solver<K>::solve(const LinearSystem<K> &ls, const Mesh &mesh, K *x) const {
		
		SolverContext ctx(ls);
	
		//DumpVector((double *)&ls.AIG.val[0], 2 * (ls.AIG.rowPtr[ls.AIG.nrow]), "AIG.val");
		//DumpVector((double *)&ls.AII.val[0], 2 * (ls.AII.rowPtr[ls.AII.nrow]), "AII.val");
		//DumpVector((double *)&ls.AGG.val[0], 2 * (ls.AIG.rowPtr[ls.AGG.nrow]), "AGG.val");

		// LU-factorize the interior-interior block of the stiffness matrix
		ctx.factorizedIStiffness = LU::Factorize(ls.AII);
		// Push the other blocks of the stiffness matrix on the GPU
		ctx.AIGDevice.reset(new CudaSparseMatrix<K>(ls.AIG));
		ctx.AGGDevice.reset(new CudaSparseMatrix<K>(ls.AGG));

		// Compute the right hand side of the Schur complement system in parallel
		std::vector<K> bt;
		computeBTilde(ctx, &bt);

		using std::tr1::placeholders::_1;
		using std::tr1::placeholders::_2;
		using std::tr1::placeholders::_3;
		MatrixProductFunctor_t schurComplement(std::tr1::bind(&Solver::multiplyBySchurComplement, this, _1, _2, _3));
		
		K *trace = x + mesh.NbLocalInterior();
		conjugateGradient(ctx, schurComplement, bt, trace, 1e-10, 100000);

		solveLocalInterior(ctx, x);
	}

	// bti = bG - AGI * AII^{-1} * bI
	template<typename K>
	void Solver<K>::computeBTilde(const SolverContext &ctx, std::vector<K> *bt) const {

		const size_t nbInterior = ctx.AIGDevice->nrow;
		const size_t nbInterface = ctx.AIGDevice->ncol;

		// xi = AII^{-1} * bI
		std::vector<K> x(nbInterior);
		ctx.factorizedIStiffness->solve(ctx.ls.loadVector, &(x[0]));

		// bti = bG - AGI * xi
		thrust::device_ptr<K> dx = ctx.workGpuBuff1;
		thrust::device_ptr<K> db = ctx.workGpuBuff2;
		thrust::copy_n(x.begin(), nbInterior, dx);
		thrust::copy_n(ctx.ls.loadVector.begin() + nbInterior, nbInterface, db);

		ctx.AIGDevice->mv(ctx.cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE, -1.0, dx.get(), 1.0, db.get());

		bt->resize(nbInterface, 0.0); // <===== Nein !
		thrust::copy_n(db, nbInterface, bt->begin());

		// bt = Sum bti
		MPI_Allreduce(MPI_IN_PLACE, &(*bt)[0], 2 * nbInterface, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	}

	// Ap = AGG * p - AGI * AII^{-1} *(AIG * p)
	template<typename K>
	void Solver<K>::multiplyBySchurComplement(const SolverContext &ctx, const std::vector<K> &p, std::vector<K> *Ap) const{

		const size_t nbInterior = ctx.AIGDevice->nrow;
		const size_t nbInterface = ctx.AIGDevice->ncol;

		thrust::device_ptr<K> dp = ctx.workGpuBuff1;	
		thrust::device_ptr<K> dx = ctx.workGpuBuff2;
		thrust::device_ptr<K> dAp = ctx.workGpuBuff3;

		// x = AIG * p
		thrust::copy_n(p.begin(), nbInterface, dp);
		ctx.AIGDevice->mv(ctx.cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, 1.0, dp.get(), 0.0, dx.get());

		// x = AII^{-1} * x
		std::vector<K> x(nbInterior);
		thrust::copy_n(dx, nbInterior, x.begin());
		ctx.factorizedIStiffness->solve(&x[0]);
		thrust::copy_n(x.begin(), nbInterior, dx);

		// Api =  AGI * x
	    ctx.AIGDevice->mv(ctx.cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE, 1.0, dx.get(), 0.0, dAp.get());

	    // Api = AGG * p - Api
	 	ctx.AGGDevice->mv(ctx.cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, 1.0, dp.get(), -1.0, dAp.get());

		thrust::copy_n(dAp, nbInterface, Ap->begin());

		// Ap = Sum Api
		MPI_Allreduce(MPI_IN_PLACE, &(*Ap)[0], 2 * nbInterface, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		//MPI_Allreduce(MPI_IN_PLACE, dAp.get(), 2 * ctx.ls.AGG.nrow, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		//thrust::copy_n(dAp, ctx.ls.AGG.nrow, Ap->begin());
	}


	// sol = AII^{-1} (b - AIG * trace)
	template<typename K>
	void Solver<K>::solveLocalInterior(const SolverContext &ctx, K *x) const {

		const size_t nbInterior = ctx.AIGDevice->nrow;
		const size_t nbInterface = ctx.AIGDevice->ncol;

		K *interior = x;
		K *trace = x + nbInterior;

		// b = b - AIG * trace
		thrust::device_ptr<K> db = ctx.workGpuBuff1;
		thrust::device_ptr<K> dtrc = ctx.workGpuBuff2;

		thrust::copy_n(ctx.ls.loadVector.begin(), nbInterior, db);
		thrust::copy_n(trace, nbInterface, dtrc);

		ctx.AIGDevice->mv(ctx.cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, -1.0, dtrc.get(), 1.0, db.get());

		// sol = AII^{-1} * b	
		thrust::copy_n(db, nbInterior, interior);
		ctx.factorizedIStiffness->solve(interior);
	}

	template<typename K>
	void Solver<K>::consolidateLocalSolutions(const Mesh &mesh, const K *localSol, int rank) const
	{
		int mpiRank, mpiSize;
		MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
		MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

		int nbInterior = mesh.NbLocalInterior();
		int totalNbInterior = mesh.NbGlobalInterior();

		// Prepare an offset array used to concatenate the solution vectors
		std::vector<int> allNbInterior(mpiRank == rank ? mpiSize : 0);
 		MPI_Gather(&nbInterior, 1, MPI_INT, &allNbInterior[0], 1, MPI_INT, rank, MPI_COMM_WORLD);

 		std::vector<int> allOffsets(mpiSize + 1);
 		allOffsets[0] = 0;
 		std::partial_sum(allNbInterior.begin(), allNbInterior.end(), &allOffsets[1]);
 		
 		// Gather all interior vertex solutions
 		std::vector<int> allVids(mpiRank == rank ? totalNbInterior : 0);
 		MPI_Gatherv(const_cast<int *>(&mesh._localToGlobal[0]), nbInterior, MPI_INT,  &allVids[0], &allNbInterior[0], &allOffsets[0], MPI_INT, rank, MPI_COMM_WORLD);

		// Gather all interior vertex global indices
		if (mpiRank == rank) {
	 		for (int i = 0; i < mpiSize; ++i) {
	 			allNbInterior[i] *= 2;
	 			allOffsets[i] *= 2;
	 		}
	 		allOffsets[mpiSize] *= 2;
	 	}

 		std::vector<K> globalISol(mpiRank == rank ? totalNbInterior : 0);
 		MPI_Gatherv(const_cast<double *>(reinterpret_cast<const double *>(localSol)), 2 * nbInterior, MPI_DOUBLE,  &globalISol[0], &allNbInterior[0], &allOffsets[0], MPI_DOUBLE, rank, MPI_COMM_WORLD);

 		// Write the global solution
 		if (mpiRank == rank) {
 			std::vector<K> solution(mesh._nbTotalVertices + 1);
 			for (int i  = 0; i < allVids.size(); ++i) {
 				solution[allVids[i]] = globalISol[i];
 			}

 			for (int i = nbInterior; i < nbInterior + mesh._nbInterfaceVertices; ++i) {
 				solution[mesh._localToGlobal[i]] = localSol[i];
 			}

			FILE *fd = fopen("lyttelton_osm_small.sol", "w");
			fprintf(fd, "MeshVersionFormatted 1\n\nDimension 2\n\nSolAtVertices\n%d\n1 1\n\n", solution.size() - 1);
			for (int i = 1; i < solution.size(); ++i) {
				fprintf(fd, "%0.12f\n", solution[i].real());
			}
			fclose(fd);
 		}
 	}

}
