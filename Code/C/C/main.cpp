#include <stdio.h>
#include <Windows.h>
#include <math.h>

#define M_2_PI 6.28318530718

typedef struct
{
    double re;
    double im;
} complex;

unsigned int reverse(register unsigned int x)
{
    x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
    x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
    x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
    x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
    return((x >> 16) | (x << 16));
}

void naive_dft(complex *x, complex *X, int N);
void fft(complex *x, complex *X, int N);
void create_W(unsigned int N);
void printTime(LARGE_INTEGER tStart, LARGE_INTEGER tStop, LARGE_INTEGER freq);
void printResult(complex *c, int n, char *str, int verified);
int verify_impulse(complex *c, int size);

int main()
{
    LARGE_INTEGER freq, tStart, tStop;
    const int N = 4096;

    create_W(4);
    getchar();
    return 0;

    /* Get ticks per second */
    QueryPerformanceFrequency(&freq);

    /* Prep data */
    complex impulse[N];
    complex res[N];
    for (int i = 0; i < N; ++i)
    {
        impulse[i] = { 0.0, 0.0 };
        res[i] = { 0.0, 0.0 };
    }
    impulse[1].re = 1.0;
    
    /* Warm up */
    printf("Warm up...\n");
    naive_dft(impulse, res, N);

    /* Run FFT */
    printf("Running...\n");
    QueryPerformanceCounter(&tStart);
    naive_dft(impulse, res, N);
    QueryPerformanceCounter(&tStop);

    /* Print results */
    printTime(tStart, tStop, freq);
    printResult(res, N, "impulse", verify_impulse(res, N));

    getchar();
    return 0;
}

/* Naive Discrete Fourier Transform, essentially as per definition */
void naive_dft(complex *x, complex *X, int N)
{
    complex sum = { 0.0, 0.0 };
    complex y = { 0.0, 0.0 };
    complex tmp = { 0.0, 0.0 };
    double theta = 1.0;
    double c1 = -M_2_PI / N;
    double c2 = 1.0;
    for (int k = 0; k < N; ++k)
    {
        sum = { 0.0, 0.0 };
        c2 = c1 * k;
        for (int n = 0; n < N; ++n)
        {
            theta = c2 * n;
            y = { cos(theta), sin(theta) };
            sum.re += x[n].re * y.re + x[n].im * y.im;
            sum.im += x[n].re * y.im + x[n].im * y.re;
        }
        X[k] = sum;
    }
}

/* Fast Fourier Transform*/
void fft(complex *c, int N)
{
    
}

void create_W(const unsigned int N)
{
    const unsigned int depth = round(log((double)N) * 1.44269504088896340736);
    unsigned int* table = (unsigned int*)malloc(sizeof(unsigned int) * N *depth);
    for (unsigned int i = 0; i < depth; ++i)
    {
        printf("\ni\tj\tshift\tlen\trev\n");
        for (unsigned int j = 0; j < N; ++j)
        {
            printf("%u\t%u\t%u\t%u\t%u\n", i, j, ((unsigned int)j >> (depth - i - 1)), (depth - i - 1), reverse((unsigned int)((unsigned int)j >> (depth - i - 1))));

            table[i * N + j] = reverse((unsigned int)((unsigned int)j >> (depth - i - 1)));
        }
    }
    /*
    for (unsigned int i = 0; i < N; ++i)
    {
        printf("%d: ", i);
        for (unsigned int j = 0; j < depth; ++j)
        {
            printf("%u ", table[j * N + i]);
        }
        printf("\n");
    }
    */
    free(table);
}

void printTime(LARGE_INTEGER tStart, LARGE_INTEGER tStop, LARGE_INTEGER freq)
{
    printf("Time (ms): %f\n", (double)(tStop.QuadPart - tStart.QuadPart) * 1000.0 / (double)freq.QuadPart);
}

void printResult(complex *c, int n, char *str, int verified)
{
    int len = n > 10 ? 10 : n;
    printf("\nResult %s:\n", str);
    for (int i = 0; i < len; ++i)
        printf("{ %f,\t%fi }\n", c[i].re, c[i].im);
    if (len != n)
        printf("...\n");
    printf("\n%s\n", verified ? "Successful" : "Error");
}

/* Quick debugging on samples */
/* The result should have constant magnitude in the transform domain. */

int verify_impulse(complex *c, int size)
{
    for (int i = 0; i < size; ++i)
        if (abs(sqrt(c[i].re * c[i].re + c[i].im * c[i].im) - 1.0) > 0.000001)
            return 0;
    return 1;
}