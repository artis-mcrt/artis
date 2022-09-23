#include "gsl_managed.h"

#include "sn3d.h"

// these functions modify their GSL implementations for CUDA compatability
// memory is allocated in managed memory instead of host memory

gsl_matrix *gsl_matrix_alloc_managed(const size_t n1, const size_t n2) {
#if CUDA_ENABLED
  gsl_block *block;
  cudaMallocManaged(&block, sizeof(gsl_block));
  block->size = sizeof(double) * n1 * n2;
  cudaMallocManaged((double **)&block->data, block->size);

  gsl_matrix *m;
  cudaMallocManaged(&m, sizeof(gsl_matrix));

  m->data = block->data;
  m->size1 = n1;
  m->size2 = n2;
  m->tda = n2;
  m->block = block;
  m->owner = 1;

  return m;
#else
  return gsl_matrix_alloc(n1, n2);
#endif
}

gsl_matrix *gsl_matrix_calloc_managed(const size_t n1, const size_t n2) {
#if CUDA_ENABLED
  gsl_matrix *mat = gsl_matrix_alloc_managed(n1, n2);
  // memset(mat->block->data, 0, sizeof(double) * n1 * n2);
  // cudamemset(mat->block->data, 0, sizeof(double) * n1 * n2);
  gsl_matrix_set_zero(mat);
  return mat;
#else
  return gsl_matrix_calloc(n1, n2);
#endif
}

__device__ __host__ double *gsl_matrix_ptr_managed(gsl_matrix *m, const size_t i, const size_t j) {
  return (double *)(m->data + (i * m->tda + j));
}

void gsl_matrix_free_managed(gsl_matrix *m) {
#if CUDA_ENABLED
  cudaFree(m->block->data);
  cudaFree(m->block);
  cudaFree(m);
#else
  gsl_matrix_free(m);
#endif
}

gsl_vector *gsl_vector_alloc_managed(const size_t n, bool readmostly) {
#if CUDA_ENABLED
  // int myGpuId;
  // cudaGetDevice(&myGpuId);

  gsl_block *block;
  cudaMallocManaged(&block, sizeof(gsl_block));

  block->size = sizeof(double) * n;
  cudaMallocManaged(&block->data, block->size);

  gsl_vector *v;
  cudaMallocManaged(&v, sizeof(gsl_vector));

  if (readmostly) {
    cudaMemAdvise(block, sizeof(gsl_block), cudaMemAdviseSetReadMostly, myGpuId);
    cudaMemAdvise(block->data, block->size, cudaMemAdviseSetReadMostly, myGpuId);
    cudaMemAdvise(v, sizeof(gsl_vector), cudaMemAdviseSetReadMostly, myGpuId);
  }

  v->data = block->data;
  v->size = n;
  v->stride = 1;
  v->block = block;
  v->owner = 1;

  return v;
#else
  return gsl_vector_alloc(n);
#endif
}

gsl_vector *gsl_vector_calloc_managed(const size_t n, bool readmostly) {
#if CUDA_ENABLED
  gsl_vector *vec = gsl_vector_alloc_managed(n, readmostly);
  gsl_vector_set_zero(vec);
  return vec;
#else
  return gsl_vector_calloc(n);
#endif
}

__device__ __host__ double gsl_vector_get_managed(gsl_vector *v, const size_t i) { return v->data[i * v->stride]; }

void gsl_vector_free_managed(gsl_vector *vec) {
#if CUDA_ENABLED
  cudaFree(vec->block->data);
  cudaFree(vec->block);
  cudaFree(vec);
#else
  gsl_vector_free(vec);
#endif
}
