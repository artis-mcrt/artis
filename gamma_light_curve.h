#ifndef GAMMA_LIGHT_CURVE_H
#define GAMMA_LIGHT_CURVE_H

int make_gamma_light_curve(void);
int gather_gamma_light_curve(int my_rank);
int write_gamma_light_curve(void);
int add_to_lc_angle(PKT *pkt_ptr);

#endif //GAMMA_LIGHT_CURVE_H
