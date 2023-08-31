#ifndef BOUNDARY_H
#define BOUNDARY_H

enum cell_boundary {
  COORD0_MIN = 101,
  COORD0_MAX = 102,
  COORD1_MIN = 103,
  COORD1_MAX = 104,
  COORD2_MIN = 105,
  COORD2_MAX = 106,
  NONE = 107,
};

double boundary_cross(struct packet *pkt_ptr, int *snext);
void change_cell(struct packet *pkt_ptr, int snext);

#endif  // BOUNDARY_H
