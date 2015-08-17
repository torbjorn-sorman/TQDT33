#ifndef MRTEST_H
#define MRTEST_H

#include "mrmath.h"

int run_simpleTest(uint32_t N);
int run_fbtest(void(*fn)(int, my_complex*, my_complex*, uint32_t), uint32_t N);
double run_test(void(*fn)(int, my_complex*, my_complex*, uint32_t), uint32_t N);
double run_2dtest(void(*fn)(int, my_complex*, my_complex*, uint32_t), my_complex **in, uint32_t N);
int run_fft2dinvtest(void(*fn)(int, my_complex*, my_complex*, uint32_t), uint32_t N);

int run_imgTest(void(*fn)(int, my_complex*, my_complex*, uint32_t), uint32_t N);
int run_fft2dTest(void(*fn)(int, my_complex*, my_complex*, uint32_t), uint32_t N);

void fft_ref(int dir, my_complex *in, my_complex *out, uint32_t N);

#endif