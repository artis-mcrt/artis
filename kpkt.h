#ifndef KPKT_H
#define KPKT_H

enum coolingtype {
  COOLINGTYPE_FF         = 880,
  COOLINGTYPE_FB         = 881,
  COOLINGTYPE_COLLEXC    = 882,
  COOLINGTYPE_COLLION    = 883,
};

struct cellhistorycoolinglist
{
  enum coolingtype type;
  int element;
  int ion;
  int level;
  int upperlevel;
};


#include "thermalbalance.h"
#include "sn3d.h"
#include "types.h"


void setup_coolinglist(void);
__host__ __device__ void calculate_cooling_rates(int modelgridindex, struct heatingcoolingrates *heatingcoolingrates);
__host__ __device__ double do_kpkt_bb(struct packet *pkt_ptr);
__host__ __device__ double do_kpkt(struct packet *pkt_ptr, double t2, int nts);
__host__ __device__ int get_coolinglistoffset(int element, int ion);

#endif //KPKT_H
