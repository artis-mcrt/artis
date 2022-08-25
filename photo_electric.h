#ifndef PHOTO_ELECTRIC_H
#define PHOTO_ELECTRIC_H

#include "types.h"
#include "cuda.h"

__host__ __device__ double sig_photo_electric(const struct packet *pkt_ptr);
__host__ __device__ double sig_pair_prod(const struct packet *pkt_ptr);
__host__ __device__ void pair_prod(struct packet *pkt_ptr);

#endif //PHOTO_ELECTRIC_H
