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

int main()
{
    double time;
    
    time = run_2dtest(mrfft, N);
    printf("\nFFT: %fs\n\n", time);

    time = run_2dtest(fft_ref, N);
    printf("\nKISS FFT: %fs\n\n", time);

    run_imgTest(mrfft, N);
    printf("... done!\n");

    getchar();
    return 0;
}