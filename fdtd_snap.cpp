#include "fdtd_vmap.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
void snapshot2d(Grid *g) {
int mm, nn;
float dim1, dim2, temp;
char filename[100];
FILE *out;

static int frame, temporalStride, startTime,
startNodeX, endNodeX, spatialStrideX,
startNodeY, endNodeY, spatialStrideY;
static char basename[80];


temporalStride=10;
startTime=0;

startNodeX=0;
endNodeX=SizeX-2;
spatialStrideX=1;

startNodeY=0;
endNodeY=SizeY-2;
spatialStrideY=1;

strcpy(basename,"nc3");


/* get snapshot if temporal conditions met */
if (Time >= startTime && (Time - startTime) % temporalStride == 0) {
sprintf(filename, "%s.%d", basename, frame++);
out = fopen(filename, "wb");
/* write dimensions to output file --
* express dimensions as floats */
dim1 = (endNodeX - startNodeX) / spatialStrideX + 1;
dim2 = (endNodeY - startNodeY) / spatialStrideY + 1;
fwrite(&dim1, sizeof(float), 1, out);
fwrite(&dim2, sizeof(float), 1, out);
/* write remaining data */
for (nn = endNodeY; nn >= startNodeY; nn -= spatialStrideY)
for (mm = startNodeX; mm <= endNodeX; mm += spatialStrideX) {
temp = (float)Ex(mm, nn);//(Hx(mm, nn) * Hx(mm, nn) + Hy(mm, nn) * Hy(mm, nn)); // store data as a floatJz(mm, nn);//
fwrite(&temp, sizeof(float), 1, out); // write the float
}
fclose(out); // close file
}
return;
}