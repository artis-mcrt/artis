#ifndef INPUT_H
#define INPUT_H

char compositionfile_hash[33];
char adatafile_hash[33];
char phixsfile_hash[33];

void input(int rank);
void read_parameterfile(int rank);
void update_parameterfile(int nts);

#endif //INPUT_H
