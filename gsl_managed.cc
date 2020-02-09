#include "gsl_managed.h"


gsl_matrix *gsl_matrix_alloc_managed(const size_t n1, const size_t n2)
{
  gsl_block *block = (gsl_block *) malloc(sizeof(gsl_block));
  block->size = sizeof(double) * n1 * n2;
  cudaMallocManaged((double **) &block->data, sizeof(double) * n1 * n2);

  gsl_matrix *mat = gsl_matrix_alloc_from_block(block, 0, n1, n2, n2);
  // mat->owner = 1;
  return mat;
}


gsl_matrix *gsl_matrix_calloc_managed(const size_t n1, const size_t n2)
{
  gsl_matrix *mat = gsl_matrix_alloc_managed(n1, n2);
  // memset(mat->block->data, 0, sizeof(double) * n1 * n2);
  // cudamemset(mat->block->data, 0, sizeof(double) * n1 * n2);
  gsl_matrix_set_zero(mat);
  return mat;
}

__device__ __host__
double *gsl_matrix_ptr_managed(gsl_matrix * m, const size_t i, const size_t j)
{
  return (double *) (m->data + (i * m->tda + j)) ;
}


void gsl_matrix_free_managed(gsl_matrix *m)
{
  if (!m->owner)
  {
    // printf("gsl_matrix_free_managed cudaFree(m->block->data)");
    cudaFree(m->block->data);
    // printf("gsl_matrix_free_managed free(m->block)");
    free(m->block);
  }
  // printf("gsl_matrix_free_managed gsl_matrix_free(m)");
  gsl_matrix_free(m);
}


gsl_vector *gsl_vector_alloc_managed(const size_t n)
{
  gsl_block *block = (gsl_block *) malloc(sizeof(gsl_block));
  block->size = sizeof(double) * n;
  cudaMallocManaged((double **) &block->data, sizeof(double) * n);

  gsl_vector *vec = gsl_vector_alloc_from_block(block, 0, n, 1);
  // vec->owner = 1;
  return vec;
}


gsl_vector *gsl_vector_calloc_managed(const size_t n)
{
  gsl_vector *vec = gsl_vector_alloc_managed(n);
  gsl_vector_set_zero(vec);
  return vec;
}


__device__ __host__
double gsl_vector_get_managed(gsl_vector *v, const size_t i)
{
  return v->data[i * v->stride];
}


void gsl_vector_free_managed(gsl_vector *vec)
{
  if (!vec->owner)
  {
    // printf("gsl_vector_free_managed cudaFree(vec->block->data);");
    cudaFree(vec->block->data);
    // printf("gsl_vector_free_managed free(vec->block);");
    free(vec->block);
  }
  // printf("gsl_vector_free_managed gsl_vector_free(vec);");
  gsl_vector_free(vec);
}
