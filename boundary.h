#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "types.h"

double boundary_cross(PKT *pkt_ptr, double tstart, int *snext);
void change_cell(PKT *pkt_ptr, int snext, bool *end_packet, double t_current, int timestep);
void change_cell_vpkt(PKT *pkt_ptr, int snext, bool *end_packet, double t_current);

#endif //BOUNDARY_H
