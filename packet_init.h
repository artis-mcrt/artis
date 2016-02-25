#ifndef PACKET_INIT_H
  #define PACKET_INIT_H

  int packet_init(int middle_iteration, int my_rank);
  int setup_packets(int pktnumberoffset);
  double fni(CELL *grid_ptr);
  double f52fe(CELL *grid_ptr);
  double f48cr(CELL *grid_ptr);
  int place_pellet(struct grid *grid_ptr, double e0, int m, int n, int pktnumberoffset);
  void write_packets(FILE *packets_file);
  void read_packets(FILE *packets_file);

#endif //PACKET_INIT_H
