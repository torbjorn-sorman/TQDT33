#ifndef MRFFT_H
#define MRFFT_H

#include "mrmath.h"

#define FORWARD_FFT 1
#define INVERSE_FFT -1

void mrfft(int dir, my_complex *x, my_complex *X, uint32_t N);
void mrfft2d(int dir, void(*fn)(int, my_complex*, my_complex*, uint32_t), my_complex **seq2d, uint32_t N);
void mrdft_naive(my_complex *x, my_complex *X, uint32_t N);

#endif