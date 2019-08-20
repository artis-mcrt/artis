#include <omp.h>
#include <mpi.h>

#ifdef _OPENMP
#define OMP_FOR_DYNAMIC _Pragma("omp for schedule(dynamic) nowait")
#define OMP_PARALLEL(x)                                             \
  int num_threads = omp_get_max_threads();                          \
  double *timer, *start_ms, *end_ms;                                \
  timer = malloc(num_threads*sizeof(double));                       \
  if (timer == NULL) abort();                                       \
  start_ms = malloc(num_threads*sizeof(double));                    \
  if (start_ms == NULL) abort();                                    \
  end_ms = malloc(num_threads*sizeof(double));                      \
  if (end_ms == NULL) abort();                                      \
                                                                    \
  _Pragma("omp parallel") {                                         \
  int tid = omp_get_thread_num();                                   \
  double start = omp_get_wtime();                                   \
  x                                                                 \
  double end = omp_get_wtime();                                     \
  timer[tid] = (end - start) * 1000.;                               \
  start_ms[tid] = start * 1000.;                                    \
  end_ms[tid] = end * 1000.;                                        \
  }                                                                 \
                                                                    \
  int rank;                                                         \
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);                              \
  for (int i = 0; i < num_threads; i++) printout("thread_timing %d  \
    %d %.12e %.12e %.12e %s\n", rank, i, timer[i], start_ms[i],     \
    end_ms[i], __FUNCTION__);                                       \
  free(timer);                                                      \
  free(start_ms);                                                   \
  free(end_ms);
#else
#define OMP_PARALLEL(x) x
#define OMP_FOR_DYNAMIC
#define OMP_FOR_DYNAMIC_REDUCTION(x) 
#endif
