#include "fdtd_vmap.h"
#include "part_init.h"
#include "part_func.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
void part_snap(Particle ph, Grid *g) {
int mm, nn, ii;
float dim1, dim2, temp;
char filename[100];
FILE *out;

static int frame, temporalStride, startTime,
startNodeX, endNodeX, spatialStrideX,
startNodeY, endNodeY, spatialStrideY;
static char basename[80];


temporalStride=10;
startTime=0;


strcpy(basename,"particle");

/* get snapshot if temporal conditions met */
if (Time >= startTime && (Time - startTime) % temporalStride == 0) {
sprintf(filename, "%s.%d", basename, frame++);
out = fopen(filename, "wb");

ii=0;
while(ph) {

if(one()*500<1.0){	
	
temp = (float)ph->x; // store data as a float

fwrite(&temp, sizeof(float), 1, out); // write the float

temp = (float)ph->y; // store data as a float

fwrite(&temp, sizeof(float), 1, out); // write the float

temp = (float)(ph->ID*1.0); // store data as a float

fwrite(&temp, sizeof(float), 1, out); // write the float

ph=ph->pNext;
ii++;
}else {ph=ph->pNext;
	   ii++;}
}
fclose(out); // close file
}
return;
}