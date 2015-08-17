#include <stdio.h>
#include <Windows.h>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <string.h>

#include "kiss_fft.h"

#include "mrmath.h"
#include "mrimage.h"
#include "mrfft.h"
#include "mrtest.h"

#define N 512 // 8192, 65536, 1048576, 2097152, 4194304, 8388608, 16777216

typedef kiss_fft_cpx my_complex;

short FFT(short int dir, long m, my_complex *);

void printTime(LARGE_INTEGER tStart, LARGE_INTEGER tStop, LARGE_INTEGER freq);
void printResult(my_complex *c, int n, char *str, int verified);

int main()
{
    run_imgTest("lena_std.ppm", mrfft, N);
    getchar();
    return 0;
}

void printTime(LARGE_INTEGER tStart, LARGE_INTEGER tStop, LARGE_INTEGER freq)
{
    printf("Time (ms): %f\n", (float)(tStop.QuadPart - tStart.QuadPart) * 1000.0 / (float)freq.QuadPart);
}

void printResult(my_complex *c, int n, char *str, int verified)
{
    int len = n > 16 ? 16 : n;
    printf("\nResult %s:\n", str);
    for (int i = 0; i < len; ++i)
        printf("{ %f,\t%fi }\n", c[i].r, c[i].i);
    if (len != n)
        printf("...\n");
}
/*
my_complex *in = (my_complex *)malloc(sizeof(my_complex) * N);
my_complex *out = (my_complex *)malloc(sizeof(my_complex) * N);

for (int i = 0; i < N; ++i)
{
in[i].r = sin(i * (M_PI / N));
in[i].i = 0.f;
}
printResult(in, N, "in", 1);
fft(1, in, out);
printResult(out, N, "in", 1);
fft(-1, out, in);
printResult(in, N, "in", 1);
*/

//my_complex **file_in = (my_complex **)malloc(sizeof(my_complex) * N * N);
//my_complex *in = (my_complex *)malloc(sizeof(my_complex) * N);
//my_complex *out = (my_complex *)malloc(sizeof(my_complex) * N);
//my_complex *out_ref = (my_complex *)malloc(sizeof(my_complex) * N);
//my_complex *out_ref2 = (my_complex *)malloc(sizeof(my_complex) * N);
/*
for (int i = 0; i < N; ++i)
{
in[i].r = in[i].i = 0.f;
}
in[1].r = 1.f;
*/
//printf("\nRunning FFT 2D...\n");
//printf("Time: %f\t\n", time = run_2dtest(1, fft, file_in));
//printf("Time: %f\t", time = run_test(fft, in, out));

//printf("\nRunning REF FFT...\n");
//printf("Time: %f\t", time_ref2 = run_test(fft_ref2, in, out_ref2));

//printf("\nRunning KISS FFT...\n");
//printf("Time: %f\t", time_ref = run_test(fft_ref, in, out_ref));

/*
printResult(out, N, "in", 1);
printResult(out_ref, N, "in", 1);
printResult(out_ref2, N, "in", 1);
*/
//compareComplex(out, out_ref2);
//compareComplex(out, out_ref);
//compareComplex(out, out_ref);
//printf("Quota: %f\t(lower is faster)\n", (time / time_ref));
//printf("Quota: %f\t(lower is faster)\n", (time / time_ref2));

//free(in);
//free(out);
