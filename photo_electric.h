#ifndef PHOTO_ELECTRIC_H
#define PHOTO_ELECTRIC_H

double sig_photo_electric(const struct packet *pkt_ptr);
double sig_pair_prod(const struct packet *pkt_ptr);
void pair_prod(struct packet *pkt_ptr);

#endif  // PHOTO_ELECTRIC_H
