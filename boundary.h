#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "types.h"

double boundary_cross(PKT *pkt_ptr, double tstart, int *snext);
int change_cell(PKT *pkt_ptr, int snext, int *end_packet, double t_current);
int change_cell_vpkt(PKT *pkt_ptr, int snext, int *end_packet, double t_current);
int locate(const PKT *pkt_ptr, double t_current);

#endif //BOUNDARY_H
