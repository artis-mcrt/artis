#ifndef BOUNDARY_H
#define BOUNDARY_H

enum cell_boundary {
  NEG_X = 101,
  POS_X = 102,
  NEG_Y = 103,
  POS_Y = 104,
  NEG_Z = 105,
  POS_Z = 106,
  NONE = 107,
};

double boundary_cross(struct packet *pkt_ptr, int *snext);
void change_cell(struct packet *pkt_ptr, int snext);

#endif  // BOUNDARY_H
