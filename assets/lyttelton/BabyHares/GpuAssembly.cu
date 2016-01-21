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
#include "GpuAssembly.h"
#include "Mesh.h"
#include "Wave.h"

#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <vector>

#include <cuComplex.h>

namespace BabyHares {
  
  namespace {
    
    extern __shared__ float shmem[];

    __host__ __device__ static __inline__ cuFloatComplex cuCsqrt(cuFloatComplex x)
    {
      float radius = cuCabsf(x);
      float cosA = x.x / radius;
      cuFloatComplex out;
      out.x = sqrt(radius * (cosA + 1.0) / 2.0);
      out.y = sqrt(radius * (1.0 - cosA) / 2.0);
      // signbit should be false if x.y is negative
      if (signbit(x.y))
        out.y *= -1.0;
      return out;
    }

    __global__
    void GatherVertexCoordinates( const int *triangles,
                                  size_t nbTriangles,
                                  int vi,
                                  const float *vertexCoords,
                                  float *gatheredCoords) {
      int myCoordId = blockIdx.x * blockDim.x + threadIdx.x, myTriangleId = myCoordId / 2;
      if (myTriangleId < nbTriangles) {
        int vertexId = triangles[3 * myTriangleId + vi];
        gatheredCoords[myCoordId] = vertexCoords[2 * vertexId + (threadIdx.x & 1)];
      }
    }

    __global__
    void Substract( const float *a, 
                    const float *b,
                    float *r, 
                    size_t length) {
      int myId = blockIdx.x * blockDim.x + threadIdx.x; 
      if (myId < length) {
        r[myId] = a[myId] - b[myId];
      }
    }
  
    __global__
    void ComputeN0s(const float *depth, float k0, float *n0s, int length) {
      int myId = blockIdx.x * blockDim.x + threadIdx.x; 
      if (myId < length) {
        float kh = k0 * depth[myId];
        float tanhkh = tanh(kh);
        n0s[myId] = 0.5 * (1 + kh * (1 - tanhkh * tanhkh) / tanhkh);
      }
    }

    __global__
    void ComputeWfs(const int *triangles, const float *depth, const cuFloatComplex *eta0, float k0, float cf, float *Wfs, int length) {
      // Wf/omega = (8 * cf * mean(abs(eta0(triangles)), 2)) ./ (3 * pi * sinh(k0 * hs).^3);
      int myId = blockIdx.x * blockDim.x + threadIdx.x; 
      if (myId < length) {
        float s = sin(k0 * depth[myId]);
        float d = 3 * M_PI * s * s * s;
        float n = 
          8.0 * cf * (
            cuCabsf(eta0[triangles[3 * myId]]) 
            + cuCabsf(eta0[triangles[3 * myId + 1]]) 
            + cuCabsf(eta0[triangles[3 * myId + 2]])
          ) / 3.;
        Wfs[myId] = n / d;
      }
    }

    __global__
    void ComputeModifiedWaveNumber(const float *Wfs, const float *n0s, float k0, cuFloatComplex *ps, int length) {
      // The modified wave number
      // ps = i * k0 * sqrt(1 - (i * Wf) ./ n0s);
      int myId = blockIdx.x * blockDim.x + threadIdx.x; 
      if (myId < length) {
        ps[myId] = cuCmulf(
                      make_cuComplex(k0, 0),
                      cuCsqrt(
                        make_cuComplex(1, -Wfs[myId] / n0s[myId])
                      )
                  );
      }
    }

    __global__
    void ComputeAreas(const float *u, const float *v, float *areas, size_t nbTriangles) {

      float *su = &shmem[0], *sv = &su[blockDim.x];

      int myCoordId = blockIdx.x * blockDim.x + threadIdx.x;
      if (myCoordId < 2 * nbTriangles) {
         su[threadIdx.x] = u[myCoordId];
         sv[threadIdx.x] = v[myCoordId];
      }

      __syncthreads();

       int myTriangleId = blockIdx.x * (blockDim.x / 2) + threadIdx.x;
       if (2 * threadIdx.x < blockDim.x && myTriangleId < nbTriangles) {
         areas[myTriangleId] = 0.5f * (su[2 * threadIdx.x] * sv[2 * threadIdx.x + 1] - su[2 * threadIdx.x + 1] * sv[2 * threadIdx.x]);
       }
    }

    __global__
    void ComputeInteriorCoeffs( const float *u,
                                const float *v, 
                                const float *w,
                                const float *areas,
                                const float *n0s,
                                const float *Wfs,
                                float k0,
                                cuFloatComplex *stiffnessOnGpu,
                                size_t nbTriangles) {

      int trianglesPerBlock = blockDim.x / 2;

      float *su = &shmem[0], *sv = &su[blockDim.x], *sw = &sv[blockDim.x];
      cuFloatComplex *sstiffness = (cuFloatComplex *)&sw[blockDim.x];

      int myCoordId = blockIdx.x * blockDim.x + threadIdx.x;
      if (myCoordId < 2 * nbTriangles) {
         su[threadIdx.x] = u[myCoordId];
         sv[threadIdx.x] = v[myCoordId];
         sw[threadIdx.x] = w[myCoordId];
      }

     __syncthreads();

     int myTriangleId = blockIdx.x * trianglesPerBlock + threadIdx.x;
     if (2 * threadIdx.x < blockDim.x && myTriangleId < nbTriangles) {

        cuFloatComplex *myStiffness = sstiffness + 9 * threadIdx.x;

        float u1 = su[2 * threadIdx.x];
        float u2 = su[2 * threadIdx.x + 1];
        float v1 = sv[2 * threadIdx.x];
        float v2 = sv[2 * threadIdx.x + 1];
        float w1 = sw[2 * threadIdx.x];
        float w2 = sw[2 * threadIdx.x + 1];

        float n0 = n0s[myTriangleId], area = areas[myTriangleId];

        float Wf = Wfs[myTriangleId];

        float prefactor = -n0 /  (4 * k0 * k0 * area);

        cuFloatComplex uu = make_cuComplex(prefactor * (u1 * u1 + u2 * u2), 0);
        cuFloatComplex uv = make_cuComplex(prefactor * (u1 * v1 + u2 * v2), 0);
        cuFloatComplex uw = make_cuComplex(prefactor * (u1 * w1 + u2 * w2), 0);
        cuFloatComplex vv = make_cuComplex(prefactor * (v1 * v1 + v2 * v2), 0);
        cuFloatComplex vw = make_cuComplex(prefactor * (v1 * w1 + v2 * w2), 0);
        cuFloatComplex ww = make_cuComplex(prefactor * (w1 * w1 + w2 * w2), 0);

        const cuFloatComplex imaginaryUnit = make_cuComplex(0, 1);

        // massDiag = (n0s - i * (Wf ./ omega)) .* areas / 6;
        cuFloatComplex mDiag = 
          cuCmulf(
            make_cuComplex(n0, -Wf),
            make_cuComplex(area / 6, 0)
          );

        cuFloatComplex m = cuCmulf(mDiag, make_cuComplex(0.5, 0));

        myStiffness[0] = cuCaddf(uu, mDiag);
        myStiffness[1] = cuCaddf(uv, m);
        myStiffness[2] = cuCaddf(uw, m);
        myStiffness[3] = cuCaddf(uv, m);
        myStiffness[4] = cuCaddf(vv, mDiag);
        myStiffness[5] = cuCaddf(vw, m);
        myStiffness[6] = cuCaddf(uw, m);
        myStiffness[7] = cuCaddf(vw, m);
        myStiffness[8] = cuCaddf(ww, mDiag);
       }

      __syncthreads();

      int coeffsPerBlock = 9 * trianglesPerBlock, totalCoeff = 9 * nbTriangles;

      for ( int myStiffnessId = blockIdx.x * coeffsPerBlock + threadIdx.x, localStiffnessId = threadIdx.x;
            localStiffnessId < coeffsPerBlock  && myStiffnessId < totalCoeff;
            myStiffnessId += blockDim.x, localStiffnessId += blockDim.x) {
            stiffnessOnGpu[myStiffnessId] = sstiffness[localStiffnessId];
      }

    }

    template <typename T>
    typename T::value_type *CudaRawPtr(T &container) {
      return  thrust::raw_pointer_cast(container.data());
    }

  }

  void ComputeInteriorCoefficientsOnGpu(const Wave &wave, 
                                        const Mesh &mesh,
                                        double cf,
                                        const std::complex<double> *eta0,
                                        std::vector<std::complex<double> > &coeffs, 
                                        std::vector<double> &n0s,
                                        std::vector<std::complex<double> > &ps) { 
    
      const size_t nbVertices = mesh._vertices.size();

      std::vector<double> vertexCoords(2 * nbVertices);
      for (size_t vi = 0; vi < nbVertices; ++vi) {
        vertexCoords[2 * vi] = mesh._vertices[vi].x;
        vertexCoords[2 * vi + 1] = mesh._vertices[vi].y;
      }

      int nbTriangles =  mesh._triangles.size();

      thrust::device_vector<float> vertexCoordsOnGpu(vertexCoords.begin(), vertexCoords.end());
      thrust::device_vector<int> triangleVidsOnGpu(3 * nbTriangles);

      cudaMemcpy2D( CudaRawPtr(triangleVidsOnGpu),
                    3 * sizeof(int),
                    &mesh._triangles[0] + offsetof(Triangle, v),
                    sizeof(Triangle),
                    3 * sizeof(int),
                    nbTriangles,
                    cudaMemcpyHostToDevice);
    


      int trianglesPerBlock = 512,  coordinatesPerBlock = 2 * trianglesPerBlock, nbBlock = 1 + nbTriangles / trianglesPerBlock;

      thrust::device_vector<float> q1(2 * nbTriangles),  q2(2 * nbTriangles),  q3(2 * nbTriangles);

      GatherVertexCoordinates<<<nbBlock, coordinatesPerBlock>>>(CudaRawPtr(triangleVidsOnGpu),
                                                                  nbTriangles,
                                                                  0, 
                                                                  CudaRawPtr(vertexCoordsOnGpu),
                                                                  CudaRawPtr(q1));
      cudaDeviceSynchronize();

      GatherVertexCoordinates<<<nbBlock, coordinatesPerBlock>>>(CudaRawPtr(triangleVidsOnGpu),
                                                                  nbTriangles,
                                                                  1, 
                                                                  CudaRawPtr(vertexCoordsOnGpu),
                                                                  CudaRawPtr(q2));
      cudaDeviceSynchronize();

      GatherVertexCoordinates<<<nbBlock, coordinatesPerBlock>>>(CudaRawPtr(triangleVidsOnGpu),
                                                                  nbTriangles,
                                                                  2, 
                                                                  CudaRawPtr(vertexCoordsOnGpu),
                                                                  CudaRawPtr(q3));
      cudaDeviceSynchronize();

     
      thrust::device_vector<float> u(2 * nbTriangles),  v(2 * nbTriangles),  w(2 * nbTriangles);
    
      Substract<<<nbBlock, 2 * trianglesPerBlock>>>(CudaRawPtr(q2), CudaRawPtr(q3), CudaRawPtr(u), u.size());
      cudaDeviceSynchronize();

      Substract<<<nbBlock, 2 * trianglesPerBlock>>>(CudaRawPtr(q3), CudaRawPtr(q1), CudaRawPtr(v), v.size());
      cudaDeviceSynchronize();

      Substract<<<nbBlock, 2 * trianglesPerBlock>>>(CudaRawPtr(q1), CudaRawPtr(q2), CudaRawPtr(w), w.size());
      cudaDeviceSynchronize();

      double k0 =  wave.waveNumber();

      thrust::device_vector<float> depthOnGpu(mesh._depth.begin(), mesh._depth.end());
      thrust::device_vector<float> n0sOnGpu(nbTriangles);
      ComputeN0s<<<nbBlock, trianglesPerBlock>>>( CudaRawPtr(depthOnGpu),
                                                  k0,
                                                  CudaRawPtr(n0sOnGpu),
                                                  nbTriangles);

      cudaDeviceSynchronize();

      thrust::device_vector<float> areasOnGpu(nbTriangles);
      ComputeAreas<<<nbBlock, coordinatesPerBlock, (2 + 2) * trianglesPerBlock * sizeof(float)>>>(CudaRawPtr(u),
                                                                                                  CudaRawPtr(v),
                                                                                                  CudaRawPtr(areasOnGpu),                                                                                                    
                                                                                                  nbTriangles);
      cudaDeviceSynchronize();

      thrust::device_vector<cuFloatComplex> eta0OnGpu(nbVertices);
      thrust::device_ptr<float> cppCplxEta0OnGpu(reinterpret_cast<float *>(CudaRawPtr(eta0OnGpu)));
      thrust::copy_n(reinterpret_cast<const double *>(eta0), 2 * nbVertices, cppCplxEta0OnGpu);
      thrust::device_vector<float> WfsOnGpu(nbTriangles);

      ComputeWfs<<<nbBlock, trianglesPerBlock>>>( CudaRawPtr(triangleVidsOnGpu), 
                                                  CudaRawPtr(depthOnGpu), 
                                                  CudaRawPtr(eta0OnGpu),
                                                  k0, 
                                                  cf, 
                                                  CudaRawPtr(WfsOnGpu),
                                                  nbTriangles);
     
      cudaDeviceSynchronize();

      thrust::device_vector<cuFloatComplex> coeffsOnGpu(9 * nbTriangles);
      ComputeInteriorCoeffs<<<nbBlock, coordinatesPerBlock, trianglesPerBlock * ((2 + 2 + 2) * sizeof(float) + 9 * sizeof(cuFloatComplex))>>>(CudaRawPtr(u),
                                                                                                                    CudaRawPtr(v),
                                                                                                                    CudaRawPtr(w),
                                                                                                                    CudaRawPtr(areasOnGpu),
                                                                                                                    CudaRawPtr(n0sOnGpu),
                                                                                                                    CudaRawPtr(WfsOnGpu),
                                                                                                                    k0,
                                                                                                                    CudaRawPtr(coeffsOnGpu),
                                                                                                                    nbTriangles);
      thrust::device_ptr<float> cppCplxCoeffsOnGpu(reinterpret_cast<float *>(CudaRawPtr(coeffsOnGpu)));
      thrust::copy_n(cppCplxCoeffsOnGpu, 2 * coeffsOnGpu.size(), reinterpret_cast<double *>(&coeffs[0]));
      thrust::copy(n0sOnGpu.begin(), n0sOnGpu.end(), n0s.begin());

      cudaDeviceSynchronize();

      thrust::device_vector<cuFloatComplex> psOnGpu(nbTriangles);
      ComputeModifiedWaveNumber<<<nbBlock, trianglesPerBlock>>>(CudaRawPtr(WfsOnGpu), CudaRawPtr(n0sOnGpu), k0, CudaRawPtr(psOnGpu), nbTriangles);

      thrust::device_ptr<float> cpppsOnGpu(reinterpret_cast<float *>(CudaRawPtr(psOnGpu)));
      thrust::copy_n(cpppsOnGpu, 2 * psOnGpu.size(), reinterpret_cast<double *>(&ps[0]));
  }
}
