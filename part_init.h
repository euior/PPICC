#ifndef _PART_INIT_H_
#define _PART_INIT_H_
#include <stdlib.h>
typedef struct ptr{
double x, y, z;
double xold, yold, zold;
double px, py, pz;
double ax, ay, az;
double m, q;
struct ptr *pNext;
int ID;
}particle, *Particle;
#endif