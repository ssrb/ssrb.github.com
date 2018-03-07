// Copyright (c) 2015, Sebastien Sydney Robert Bigot
// All rights reserved.
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <vector>

#include <cuComplex.h>

#include <cstdio>

#define SYNC(EXPR) do {EXPR; cudaDeviceSynchronize();} while(false)
//#define SYNC(EXPR) EXPR

typedef int32_t color_t;

namespace {
  
  extern __shared__ float shmem[];

  template <typename T>
  class KernelArray
  {  
  public:
    KernelArray(thrust::device_vector<T> dVec)
     :  _array(thrust::raw_pointer_cast(dVec.data())),
        _size(dVec.size())
    {}
        
    __device__
    T &operator[](int i)
    {
      return _array[i];
    }

    __device__
    const T &operator[](int i) const
    {
      return _array[i];
    }

    __host__ __device__
    size_t size() const
    {
      return _size;
    }

  private:
    T* const _array;
    const size_t _size;
  };
   

  // This kernel gather the ith (i=0,1,2) vertices x and y coords for all triangles:
  __global__
  void GatherVertexCoordinates( 
    const KernelArray<int> triangles,                                                               
    const KernelArray<float> vertexCoords,
    KernelArray<float> gatheredxs,
    KernelArray<float> gatheredys,
    int vi
  ) {
    const int myTriangleId = blockIdx.x * blockDim.x + threadIdx.x;
    const int nbTriangle = triangles.size() / 3;
    if (myTriangleId < nbTriangle) {
      // Lookup the global index of the ith vertex of my triangle 
      int vertexId = triangles[3 * myTriangleId + vi];
      // Lookup the coord in the vertex coords array
      gatheredxs[myTriangleId] = vertexCoords[2 * vertexId];
      gatheredys[myTriangleId] = vertexCoords[2 * vertexId + 1];
    }
  }

  __global__
  void Assemble(
    const KernelArray<float> coeffs, 
    const KernelArray<color_t> colors, 
    const KernelArray<int> tt, 
    KernelArray<float> result,
    color_t color) {
    const int myTriangleId = blockIdx.x * blockDim.x + threadIdx.x;
    if (myTriangleId < colors.size() && colors[myTriangleId] == color) {
      const int *mytt = &tt[9 * myTriangleId];
      const float *myCoeffs = &coeffs[9 * myTriangleId];
      result[mytt[0]] += myCoeffs[0];
      result[mytt[1]] += myCoeffs[1];
      result[mytt[2]] += myCoeffs[2];
      result[mytt[3]] += myCoeffs[3];
      result[mytt[4]] += myCoeffs[4];
      result[mytt[5]] += myCoeffs[5];
      result[mytt[6]] += myCoeffs[6];
      result[mytt[7]] += myCoeffs[7];
      result[mytt[8]] += myCoeffs[8];
    }
  }

  __global__
  void ComputeElementsCoeffs( 
    const KernelArray<float> q1xs,
    const KernelArray<float> q1ys, 
    const KernelArray<float> q2xs,
    const KernelArray<float> q2ys,
    const KernelArray<float> q3xs, 
    const KernelArray<float> q3ys,
    KernelArray<float> coefficients
  ) {
 
    const size_t nbTriangles = q1xs.size();
    const int myTriangleId = blockIdx.x * blockDim.x + threadIdx.x;
    float *shcoeffs = &shmem[9 * threadIdx.x];
    if (myTriangleId < nbTriangles) {

      float ux = q2xs[myTriangleId] - q3xs[myTriangleId];
      float uy = q2ys[myTriangleId] - q3ys[myTriangleId];
      float vx = q3xs[myTriangleId] - q1xs[myTriangleId];
      float vy = q3ys[myTriangleId] - q1ys[myTriangleId];
      float wx = q1xs[myTriangleId] - q2xs[myTriangleId];
      float wy = q1ys[myTriangleId] - q2ys[myTriangleId];

      float area = 0.5f * (ux * vy - uy * vx);
      float prefactor = -1 /  (4 * area);

      float uu = prefactor * (ux * ux + uy * uy);
      float uv = prefactor * (ux * vx + uy * vy);
      float uw = prefactor * (ux * wx + uy * wy);
      float vv = prefactor * (vx * vx + vy * vy);
      float vw = prefactor * (vx * wx + vy * wy);
      float ww = prefactor * (wx * wx + wy * wy);

      float mDiag = area / 6, m = mDiag / 2;

      shcoeffs[0] = uu + mDiag;
      shcoeffs[1] = uv + m;
      shcoeffs[2] = uw + m;
      shcoeffs[3] = uv + m;
      shcoeffs[4] = vv + mDiag;
      shcoeffs[5] = vw + m;
      shcoeffs[6] = uw + m;
      shcoeffs[7] = vw + m;
      shcoeffs[8] = ww + mDiag;
    }

    __syncthreads();

    // Coalesced write
    int coeffsPerBlock = 9 * blockDim.x, totalCoeff = 9 * nbTriangles;
    for (int local = threadIdx.x, global = blockIdx.x * coeffsPerBlock + threadIdx.x;
          local < coeffsPerBlock && global < totalCoeff;
          global += blockDim.x, local += blockDim.x) {
          coefficients[global] = shmem[local];
    }

  }

}

extern "C"
{
  void do_global_assembly_on_gpu(
    size_t nv,
    size_t nt,
    const double *vs,
    const int32_t *ts,
    const color_t *colors,
    const int32_t *rowptr,
    const int32_t *colidx,
    const int32_t *tt,
    double *coeffs
  ) 
  { 
    using namespace thrust;

    // Push triangles and vertices on the GPU
    device_vector<float> vertexCoordsOnGpu(vs, &vs[2 * nv]);
    device_vector<int> triangleVidsOnGpu(ts, &ts[3 * nt]);
  
    int trianglesPerBlock = 512,  nbBlock = 1 + nt / trianglesPerBlock;

    // Pack triangle coordinates xs, ys for each edges u, v, w
    device_vector<float> q1xs(nt), q1ys(nt), q2xs(nt), q2ys(nt), q3xs(nt), q3ys(nt);
    SYNC((
      GatherVertexCoordinates<<<nbBlock, trianglesPerBlock>>>(
        triangleVidsOnGpu,                                                                                                                         
        vertexCoordsOnGpu,
        q1xs,
        q1ys,
        0
      )
    ));

    SYNC((
      GatherVertexCoordinates<<<nbBlock, trianglesPerBlock>>>(
        triangleVidsOnGpu,                                                              
        vertexCoordsOnGpu,
        q2xs,
        q2ys,
        1
      )
    ));

    SYNC((
      GatherVertexCoordinates<<<nbBlock, trianglesPerBlock>>>(
        triangleVidsOnGpu,
        vertexCoordsOnGpu,
        q3xs,
        q3ys,
        2
      )
    ));

    device_vector<float> coeffsOnGpu(9 * nt);
    size_t shmem = 9 * trianglesPerBlock * sizeof(float);
    SYNC((ComputeElementsCoeffs<<<nbBlock, trianglesPerBlock, shmem>>>(q1xs, q1ys, q2xs, q2ys, q3xs, q3ys, coeffsOnGpu)));

    device_vector<int> colorsOnGpu(colors, &colors[nt]);
    device_vector<int> ttOnGpu(tt, &tt[9 * nt]);
    device_vector<float> result(rowptr[nv]);
    fill(result.begin(), result.end(), 0.);
    device_vector<int>::iterator maxIter = max_element(colorsOnGpu.begin(), colorsOnGpu.end());
    int nbColor = 1 + *maxIter;
    for (int color = 0; color < nbColor; ++color)
    {
      SYNC((
        Assemble<<<nbBlock, trianglesPerBlock>>>( 
          coeffsOnGpu,
          colorsOnGpu,
          ttOnGpu,
          result,
          color
        )
      ));
    }

    copy(coeffsOnGpu.begin(), coeffsOnGpu.end(), coeffs);  
  }
}