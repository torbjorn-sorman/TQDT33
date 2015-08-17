#ifndef MRIMAGE_H
#define MRIMAGE_H

#include "mrmath.h"

void imgToComplex(unsigned char *img, my_complex **com, uint32_t N);
void minRange(double *m, double *r, my_complex **com, uint32_t N);
void fftToImg(my_complex **com, unsigned char *img, uint32_t N);
void fftShift(unsigned char *in, unsigned char *out, uint32_t N);
int writeppm(char *filename, int width, int height, unsigned char *data);
unsigned char *readppm(char *filename, int *width, int *height);

#endif