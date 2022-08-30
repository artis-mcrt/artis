#ifndef __GSL_MANAGED_H__
#define __GSL_MANAGED_H__

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>
#include <stddef.h>

#include "cuda.hpp"

gsl_matrix *gsl_matrix_alloc_managed(const size_t n1, const size_t n2);
gsl_matrix *gsl_matrix_calloc_managed(const size_t n1, const size_t n2);
void gsl_matrix_free_managed(gsl_matrix *m);
__host__ __device__ double *gsl_matrix_ptr_managed(gsl_matrix *m, const size_t i, const size_t j);

gsl_vector *gsl_vector_alloc_managed(const size_t n, bool readmostly);
gsl_vector *gsl_vector_calloc_managed(const size_t n, bool readmostly);
__host__ __device__ double gsl_vector_get_managed(gsl_vector *vec, const size_t i);
void gsl_vector_free_managed(gsl_vector *m);

#endif /* __GSL_MANAGED_H__ */
