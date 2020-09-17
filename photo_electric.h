#ifndef PHOTO_ELECTRIC_H
#define PHOTO_ELECTRIC_H

#include "types.h"

double sig_photo_electric(const PKT *pkt_ptr);
double sig_pair_prod(const PKT *pkt_ptr);
void pair_prod(PKT *pkt_ptr);

#endif //PHOTO_ELECTRIC_H
