#ifndef UPDATE_PACKETS_H
#define UPDATE_PACKETS_H

int update_packets(int nts);
void update_cell(int cellnumber);
int compare_packets_byposition(const void *p1, const void *p2);
int compare_packets_bymodelgridposition(const void *p1, const void *p2);
int compare_packets_bymodelgriddensity(const void *p1, const void *p2);

#endif //UPDATE_PACKETS_H
