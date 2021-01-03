#ifndef PHOTO_ELECTRIC_H
#define PHOTO_ELECTRIC_H

#include "types.h"

__host__ __device__ double sig_photo_electric(const PKT *pkt_ptr);
__host__ __device__ double sig_pair_prod(const PKT *pkt_ptr);
__host__ __device__ void pair_prod(PKT *pkt_ptr);

#endif //PHOTO_ELECTRIC_H
