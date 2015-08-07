#include <stdio.h>
#include <Windows.h>
#include <math.h>

#define M_2_PI 6.28318530718
#define N 4096

typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;

struct sc
{
    double re;
    double im;
};

struct complex
{
    double r;
    double i;

    complex& operator=(const complex& a)
    {
        r = a.r;
        i = a.i;
        return *this;
    }

    complex operator+(const complex& a) const
    {
        return{ a.r + r, a.i + i };
    }

    complex& operator+=(const complex& a)
    {
        r += a.r;
        i += a.i;
        return *this;
    }

    complex operator*(const complex& a) const
    {
        return{ a.r * r + a.i * i, a.i * r + a.r * i };
    }

    bool operator==(const complex& a) const
    {
        return a.r == r && a.i == i;
    }
};

uint32_t
reverse(uint32_t x, uint32_t l)
{
    x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
    x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
    x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
    x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
    return((x >> 16) | (x << 16)) >> (32 - l);

}

uint32_t nk_mod_n(uint32_t n, uint32_t b, uint32_t d);
complex complex_pow(complex c, uint32_t p);
complex pow_complex(complex c, uint32_t p);
const int tab32[32] = { 0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30, 8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31 };
int log2_32(uint32_t value);
void naive_dft(complex *x, complex *X);
void fft(complex *x, complex *X);
void printTime(LARGE_INTEGER tStart, LARGE_INTEGER tStop, LARGE_INTEGER freq);
void printResult(complex *c, int n, char *str, int verified);
int verify_impulse(complex *c, int size);
void compareComplex(complex *c1, complex *c2, double t1, double t2);

int main()
{
    LARGE_INTEGER freq, tStart, tStop;
    double time_FFT, time_NDFT;

    /* Get ticks per second */
    QueryPerformanceFrequency(&freq);

    /* Prep data */
    complex impulse[N];
    complex res_FFT[N];
    complex res_NDFT[N];
    for (int i = 0; i < N; ++i)
    {
        impulse[i] = { 0.0, 0.0 };
        res_FFT[i] = { 0.0, 0.0 };
        res_NDFT[i] = { 0.0, 0.0 };
    }
    impulse[1].r = 1.0;

    /* FFT */
    printf("\nWarm up FFT...\n");
    fft(impulse, res_FFT);    
    printf("Running...\n");
    QueryPerformanceCounter(&tStart);    
    fft(impulse, res_FFT);
    QueryPerformanceCounter(&tStop);
    printTime(tStart, tStop, freq);
    time_FFT = (double)(tStop.QuadPart - tStart.QuadPart) * 1000.0 / (double)freq.QuadPart;

    /* Naive DFT */
    printf("\nWarm up naive DFT...\n");
    naive_dft(impulse, res_NDFT);
    printf("Running...\n");
    QueryPerformanceCounter(&tStart);
    naive_dft(impulse, res_NDFT);
    QueryPerformanceCounter(&tStop);
    printTime(tStart, tStop, freq);
    time_NDFT = (double)(tStop.QuadPart - tStart.QuadPart) * 1000.0 / (double)freq.QuadPart;

    //printResult(res, N, "impulse", verify_impulse(res, N));
    compareComplex(res_FFT, res_NDFT, time_FFT, time_NDFT);
    getchar();
    return 0;
}

/* Naive Discrete Fourier Transform, essentially as per definition */
void naive_dft(complex *x, complex *X)
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
            sum += (x[n] * y);
        }
        x[k] = sum;
    }
}

/* Fast Fourier Transform */
/*
3.
In the l-th array, node k (in binary form k[y-1], ..., k1, k0) has a solid line drawn to it from a node in the (l - 1)th array.
The address of the node in the (l - 1)th array is the same as node k, except that bit k[y l] must be a one.

The dashed line comes from a node in the (l - 1)th array whose address is the same, but bit k[y l] must be a zero.
*/

void fft(complex *x, complex *X)
{
    complex tmp[N];
    const uint32_t depth = log2_32(N);
    const uint32_t n_half = N / 2;
    complex W = { cos(-M_2_PI / N), sin(-M_2_PI / N) };
    uint32_t bit = 0;
    complex upper, lower;
    uint32_t u, l;
    uint32_t dist, offset;    
    for (uint32_t n = 0; n < N; ++n)
        tmp[n] = x[n];
    dist = N;
    for (uint32_t k = 0; k < depth; ++k)
    {
        bit = depth - 1 - k;
        dist = (n_half >> k);
        offset = 0;
        for (uint32_t n = 0; n < n_half; ++n)
        {
            offset += (n >= (dist + offset)) * dist * 2;
            l = (n & ~(1 << bit)) + offset;
            u = (n | 1 << bit) + offset;
            lower = tmp[l] + complex_pow(W, nk_mod_n(l, bit, depth)) * tmp[u];
            upper = tmp[l] + complex_pow(W, nk_mod_n(u, bit, depth)) * tmp[u];
            tmp[l] = lower;
            tmp[u] = upper;
        }
    }
    for (int n = 0; n < N; ++n)
    {
        X[n] = tmp[reverse(n, depth)];
    }
}

void printTime(LARGE_INTEGER tStart, LARGE_INTEGER tStop, LARGE_INTEGER freq)
{
    printf("Time (ms): %f\n", (double)(tStop.QuadPart - tStart.QuadPart) * 1000.0 / (double)freq.QuadPart);
}

void printResult(complex *c, int n, char *str, int verified)
{
    int len = n > 16 ? 16 : n;
    printf("\nResult %s:\n", str);
    for (int i = 0; i < len; ++i)
        printf("{ %f,\t%fi }\n", c[i].r, c[i].i);
    if (len != n)
        printf("...\n");
    printf("\n%s\n", verified ? "Successful" : "Error");
}

/* Quick debugging on samples */
/* The result should have constant magnitude in the transform domain. */
int verify_impulse(complex *c, int size)
{
    for (int i = 0; i < size; ++i)
        if (abs(sqrt(c[i].r * c[i].r + c[i].i * c[i].i) - 1.0) > 0.000001)
            return 0;
    return 1;
}

uint32_t nk_mod_n(uint32_t n, uint32_t b, uint32_t d)
{
    uint32_t x = (n >> b);
    x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
    x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
    x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
    x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
    return((x >> 16) | (x << 16)) >> (32 - d);
}

int log2_32(uint32_t value)
{
    value |= value >> 1; value |= value >> 2; value |= value >> 4; value |= value >> 8; value |= value >> 16;
    return tab32[(uint32_t)(value * 0x07C4ACDD) >> 27];
}

complex complex_pow(complex c, uint32_t p)
{
    double theta = acos(c.r) * p;
    return{ cos(theta), sin(theta) };
}

complex pow_complex(complex c, uint32_t p)
{
    complex a = c;
    for (uint32_t i = 1; i < p; ++i)
        a = a * c;
    return a;
}

void compareComplex(complex *c1, complex *c2, double t1, double t2)
{
    int res = 0;
    for (int i = 0; i < N; ++i)
    {
        if (c1[i].r != c2[i].r || c1[i].i != c2[i].i)
        {
            res = 1;
            break;
        }
    }
    printf("\n%s\n", res ? "EQUAL" : "NOT EQUAL");
    printf("\nImproved %f times.", (t2 / t1));
    printf("\nTheoretical limit %f times.", (((double)N * N) / ((double)N * log2_32(N))));
}