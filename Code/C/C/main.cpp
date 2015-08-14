#include <stdio.h>
#include <Windows.h>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <string.h>

#include "kiss_fft.h"

#define M_2_PI 6.28318530718f
#define M_PI 3.14159265359f
#define N 512 // 8192, 65536, 1048576, 2097152, 4194304, 8388608, 16777216

typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef kiss_fft_cpx my_complex;

static const int tab32[32] =
{
    0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30, 8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
};

static const uint32_t revTbl256[] =
{
    0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
    0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
    0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
    0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
    0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
    0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
    0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
    0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
    0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
    0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
    0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
    0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
    0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
    0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
    0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
    0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
};

#define reverseBits(X,L) (((revTbl256[X & 0xff] << 24) | (revTbl256[(X >> 8) & 0xff] << 16) | (revTbl256[(X >> 16) & 0xff] << 8) | (revTbl256[(X >> 24) & 0xff])) >> L)
#define C_MUL(R,A,B) R.r = A.r * B.r + A.i * B.i; R.i = A.r * B.i + A.i * B.r
#define C_MUL_ADD(R,A,B,C) R.r = C.r + A.r * B.r + A.i * B.i; R.i = C.i + A.r * B.i + A.i * B.r
#define C_ABS(A) sqrt(A.r * A.r + A.i * A.i)
#define C_SINCOS(C,A) C.i = sin(A); C.r = cos(A)

uint32_t reverseBitsLowMem(uint32_t x, uint32_t l);
int log2_32(uint32_t value);
void naive_dft(my_complex *x, my_complex *X);

void fft(int dir, my_complex *x, my_complex *X);
void fft_ref(int dir, my_complex *in, my_complex *out);
void fft_ref2(int dir, my_complex *x, my_complex *X);
short FFT(short int dir, long m, my_complex *);
void fft2d(int dir, void(*fn)(int, my_complex*, my_complex*), my_complex **seq2d);

void printTime(LARGE_INTEGER tStart, LARGE_INTEGER tStop, LARGE_INTEGER freq);
void printResult(my_complex *c, int n, char *str, int verified);
void compareComplex(my_complex *c1, my_complex *c2);

double run_test(int dir, void(*fn)(int, my_complex*, my_complex*), my_complex *in, my_complex *out);
double run_2dtest(int dir, void(*fn)(int, my_complex*, my_complex*), my_complex **in);

int writeppm(char *filename, int width, int height, unsigned char *data);
unsigned char *readppm(char *filename, int *width, int *height);

void imgToComplex(unsigned char *img, my_complex **com);
void realToImg(my_complex **com, unsigned char *img);
void fftToImg(my_complex **com, unsigned char *img);

int main()
{

    //double time;//, time_ref;//, time_ref2;
    char *filename = "lena_std.ppm";

    int n, m;
    unsigned char *image = readppm(filename, &n, &m);

    if (n != N || m != N)
    {
        printf("Image size not square and pw of 2.\n");
        getchar();
        return 0;
    }

    unsigned char *grey = (unsigned char *)malloc(sizeof(unsigned char) * N * N * 3);
    my_complex **file_in = (my_complex **)malloc(sizeof(my_complex) * N);
    for (int i = 0; i < N; ++i)
        file_in[i] = (my_complex *)malloc(sizeof(my_complex) * N);

    imgToComplex(image, file_in);
    realToImg(file_in, grey);
    writeppm("greyscale.ppm", N, N, grey);

    fft2d(1, fft, file_in);
    fftToImg(file_in, grey);
    writeppm("magnitude.ppm", N, N, grey);

    fft2d(-1, fft, file_in);
    realToImg(file_in, grey);
    writeppm("fftToGrey.ppm", N, N, grey);

    for (int i = 0; i < N; ++i)
        free(file_in[i]);
    free(file_in);
    free(grey);
    printf("Done...\n");
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
    getchar();
    return 0;
}

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

void imgToComplex(unsigned char *img, my_complex **com)
{
    float r, g, b, intensity;
    uint32_t px;
    for (uint32_t y = 0; y < N; ++y)
    {
        for (uint32_t x = 0; x < N; ++x)
        {
            px = y * N * 3 + x * 3;
            r = img[px];
            g = img[px + 1];
            b = img[px + 2];
            intensity = ((r + g + b) / 3.f) / 255.f;
            com[y][x].r = intensity;
            com[y][x].i = 0.f;
        }
    }
}

void minRange(double *m, double *r, my_complex **com)
{
    uint32_t px;
    double mi = DBL_MAX;
    double ma = DBL_MIN;
    for (uint32_t y = 0; y < N; ++y)
    {
        for (uint32_t x = 0; x < N; ++x)
        {
            mi = min(mi, com[y][x].r);
            ma = max(ma, com[y][x].r);
        }
    }
    printf("min: %f, max: %f\n", mi, ma);
    *m = mi;
    *r = ma - mi;
}

void realToImg(my_complex **com, unsigned char *img)
{
    uint32_t px;
    double magnitude, vmin, range;
    minRange(&vmin, &range, com);
    printf("min: %f, range: %f\n", vmin, range);
    for (uint32_t y = 0; y < N; ++y)
    {
        for (uint32_t x = 0; x < N; ++x)
        {
            px = y * N * 3 + x * 3;
            img[px] = img[px + 1] = img[px + 2] = (unsigned char)(255.0 * com[y][x].r);
        }
    }
}

void fftToImg(my_complex **com, unsigned char *img)
{
    uint32_t px;
    double magnitude;
    for (uint32_t y = 0; y < N; ++y)
    {
        for (uint32_t x = 0; x < N; ++x)
        {
            px = y * N * 3 + x * 3;
            magnitude = sqrt(com[y][x].r *com[y][x].r + com[y][x].i *com[y][x].i);
            img[px] = img[px + 1] = img[px + 2] = (unsigned char)(255.0 * magnitude);
        }
    }
}

void fft2d(int dir, void(*fn)(int, my_complex*, my_complex*), my_complex **seq2d)
{
    const uint32_t depth = log2_32(N);
    uint32_t row, col;
    my_complex *seq, *out;

    /* Transform the rows */
    seq = (my_complex *)malloc(N * sizeof(my_complex));
    out = (my_complex *)malloc(N * sizeof(my_complex));

    for (row = 0; row < N; ++row)
    {
        for (col = 0; col < N; ++col)
            seq[col] = seq2d[row][col];

        fn(dir, seq, out);

        for (col = 0; col < N; ++col)
            seq2d[row][col] = out[col];
    }
    for (col = 0; col < N; ++col)
    {
        for (row = 0; row < N; ++row)
            seq[row] = seq2d[row][col];

        fn(dir, seq, out);

        for (row = 0; row < N; ++row)
            seq2d[row][col] = out[row];
    }
    free(seq);
}

void fft_ref(int dir, my_complex *in, my_complex *out)
{
    kiss_fft_cfg cfg = kiss_fft_alloc(N, (dir == -1), 0, 0);
    kiss_fft(cfg, in, out);
    free(cfg);
}

void fft_ref2(int dir, my_complex *in, my_complex *out)
{
    for (int i = 0; i < N; ++i)
    {
        out[i] = in[i];
    }
    FFT(dir, log2_32(N), out);
}

/* Naive Fast Fourier Transform */
void fft(int dir, my_complex *x, my_complex *X)
{
    const uint32_t depth = log2_32(N);
    const uint32_t n2 = (N / 2);
    my_complex *W = (my_complex *)malloc(sizeof(my_complex) * N);
    my_complex *tmp = (my_complex *)malloc(sizeof(my_complex) * N);
    my_complex tmp_u, tmp_l;
    uint32_t trail, bit, u, l, p, dist, dist_2, offset;
    float ang, w_angle;

    w_angle = dir * -(M_2_PI / N);
    trail = 32 - depth;
    bit = 0;
    for (uint32_t n = 0; n < N; ++n)
    {
        tmp[n] = x[n];
    }
    for (uint32_t n = 0; n <= N; ++n)
    {
        ang = w_angle * n;
        C_SINCOS(W[n], ang);
    }
    dist = N;
    for (uint32_t k = 0; k < depth; ++k)
    {
        bit = depth - 1 - k;
        dist_2 = dist;
        dist = dist >> 1;
        offset = 0;
        for (uint32_t n = 0; n < n2; ++n)
        {
            offset += (n >= (dist + offset)) * dist_2;
            l = (n & ~(1 << bit)) + offset;
            u = l + dist;
            tmp_l = tmp[l];
            tmp_u = tmp[u];
            p = (l >> bit);
            C_MUL_ADD(tmp[l], W[p], tmp_u, tmp_l);
            p = (u >> bit);
            C_MUL_ADD(tmp[u], W[p], tmp_u, tmp_l);
        }
    }
    if (dir == 1)
    {
        for (uint32_t n = 0; n < N; ++n)
        {
            X[n].r = tmp[n].r / N;
            X[n].i = tmp[n].i / N;
        }
    }
    //free(tmp);
    //free(W);
}

/*
This computes an in-place complex-to-complex FFT
x and y are the real and imaginary arrays of 2^m points.
dir =  1 gives forward transform
dir = -1 gives reverse transform

Fails in the long run since the error accumulates... N = 512, error < 0.00001
*/
short FFT(short int dir, long m, my_complex *x)
{
    long n, i, i1, j, k, i2, l, l1, l2;
    float c1, c2, t1, t2, u1, u2, z;
    my_complex tx;
    /* Calculate the number of points */
    n = N;
    /* Do the bit reversal */
    i2 = n >> 1;
    j = 0;
    for (i = 0; i < n - 1; i++) {
        if (i < j) {
            tx = x[i];
            x[i] = x[j];
            x[j] = tx;
        }
        k = i2;
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }
    /* Compute the FFT */
    c1 = -1.f;
    c2 = 0.f;
    l2 = 1;
    for (l = 0; l < m; l++) {
        l1 = l2;
        l2 <<= 1;
        u1 = 1.f;
        u2 = 0.f;
        for (j = 0; j < l1; j++) {
            for (i = j; i < n; i += l2) {
                i1 = i + l1;
                t1 = u1 * x[i1].r - u2 * x[i1].i;
                t2 = u1 * x[i1].i + u2 * x[i1].r;
                x[i1].r = x[i].r - t1;
                x[i1].i = x[i].i - t2;
                x[i].r += t1;
                x[i].i += t2;
            }
            z = u1 * c1 - u2 * c2;
            u2 = u1 * c2 + u2 * c1;
            u1 = z;
        }
        c2 = sqrt((1.f - c1) / 2.f);
        if (dir == 1)
            c2 = -c2;
        c1 = sqrt((1.f + c1) / 2.f);
    }
    return(TRUE);
}

double run_test(void(*fn)(int, my_complex*, my_complex*), my_complex *in, my_complex *out)
{
    LARGE_INTEGER freq, tStart, tStop;
    double m = DBL_MAX;
    QueryPerformanceFrequency(&freq);
    for (int i = 0; i < 100; ++i)
    {
        QueryPerformanceCounter(&tStart);
        fn(0, in, out);
        QueryPerformanceCounter(&tStop);
        m = min(m, (double)(tStop.QuadPart - tStart.QuadPart) * 1000.0 / (float)freq.QuadPart);
    }
    return m;
}

double run_2dtest(int dir, void(*fn)(int, my_complex*, my_complex*), my_complex **in)
{
    LARGE_INTEGER freq, tStart, tStop;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&tStart);
    fft2d(dir, fn, in);
    QueryPerformanceCounter(&tStop);
    return (double)(tStop.QuadPart - tStart.QuadPart) * 1000.0 / (float)freq.QuadPart;
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

void compareComplex(my_complex *c1, my_complex *c2)
{
    double m = 0.00001;
    double max_r = -DBL_MAX;
    double max_i = -DBL_MAX;
    for (int i = 0; i < N; ++i)
    {
        max_r = max(abs((double)c1[i].r - (double)c2[i].r), max_r);
        max_i = max(abs((double)c1[i].i - (double)c2[i].i), max_i);
    }
    if ((max_r > m || max_i > m))
        printf("\nNOT EQUAL\nDiff: (%f, %f)\n", max_r, max_i);
    else
        printf("\nEQUAL\n");
}

int log2_32(uint32_t value)
{
    value |= value >> 1; value |= value >> 2; value |= value >> 4; value |= value >> 8; value |= value >> 16;
    return tab32[(uint32_t)(value * 0x07C4ACDD) >> 27];
}

uint32_t reverseBitsLowMem(uint32_t x, uint32_t l)
{
    x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
    x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
    x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
    x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
    return((x >> 16) | (x << 16)) >> (32 - l);
}

/* Naive Discrete Fourier Transform, essentially as per definition */
void naive_dft(my_complex *x, my_complex *X)
{
    float real, img;
    my_complex y = { 0.0, 0.0 };
    float re, im;
    my_complex tmp = { 0.0, 0.0 };
    float theta = 1.0;
    float c1 = -M_2_PI / N;
    float c2 = 1.0;
    for (int k = 0; k < N; ++k)
    {
        real = 0.0;
        img = 0.0;
        c2 = c1 * k;
        for (int n = 0; n < N; ++n)
        {
            theta = c2 * n;
            re = cos(theta);
            im = sin(theta);
            real += x[n].r * re + x[n].i * im;
            img += x[n].r * im + x[n].i * re;
        }
        x[k].r = real;
        x[k].i = img;
    }
}

int writeppm(char *filename, int width, int height, unsigned char *data)
{
    FILE *fp;
    int error = 1;
    int i, h, v;

    if (filename != NULL)
    {
        fopen_s(&fp, filename, "w");

        if (fp != NULL)
        {
            // Write PPM file
            // Header	
            fprintf(fp, "P3\n");
            fprintf(fp, "# written by Ingemars PPM writer\n");
            fprintf(fp, "%d %d\n", width, height);
            fprintf(fp, "%d\n", 255); // range

            // Data
            for (v = height - 1; v >= 0; v--)
            {
                for (h = 0; h < width; h++)
                {
                    i = (width*v + h) * 3; // assumes rgb, not rgba
                    fprintf(fp, "%d %d %d ", data[i], data[i + 1], data[i + 2]);
                }
                fprintf(fp, "\n"); // range
            }

            if (fwrite("\n", sizeof(char), 1, fp) == 1)
                error = 0; // Probable success
            fclose(fp);
        }
    }
    return(error);
}

unsigned char *readppm(char *filename, int *width, int *height)
{
    FILE *fd;
    int  k;//, nm;
    char c;
    int i, j;
    char b[100];
    //float s;
    int red, green, blue;
    long numbytes;//, howmuch;
    int n;
    int m;
    unsigned char *image;

    fopen_s(&fd, filename, "rb");
    if (fd == NULL)
    {
        printf("Could not open %s\n", filename);
        return NULL;
    }
    c = getc(fd);
    if (c == 'P' || c == 'p')
        c = getc(fd);

    if (c == '3')
    {
        printf("%s is a PPM file (plain text version)\n", filename);

        // NOTE: This is not very good PPM code! Comments are not allowed
        // except immediately after the magic number.
        c = getc(fd);
        if (c == '\n' || c == '\r') // Skip any line break and comments
        {
            c = getc(fd);
            while (c == '#')
            {
                fscanf_s(fd, "%[^\n\r] ", b);
                printf("%s\n", b);
                c = getc(fd);
            }
            ungetc(c, fd);
        }
        fscanf_s(fd, "%d %d %d", &n, &m, &k);

        printf("%d rows  %d columns  max value= %d\n", n, m, k);

        numbytes = n * m * 3;
        image = (unsigned char *)malloc(numbytes);
        if (image == NULL)
        {
            printf("Memory allocation failed!\n");
            return NULL;
        }
        for (i = m - 1; i >= 0; i--) for (j = 0; j < n; j++) // Important bug fix here!
        { // i = row, j = column
            fscanf_s(fd, "%d %d %d", &red, &green, &blue);
            image[(i*n + j) * 3] = red * 255 / k;
            image[(i*n + j) * 3 + 1] = green * 255 / k;
            image[(i*n + j) * 3 + 2] = blue * 255 / k;
        }
    }
    else
        if (c == '6')
        {
            printf("%s is a PPM file (raw version)!\n", filename);

            c = getc(fd);
            if (c == '\n' || c == '\r') // Skip any line break and comments
            {
                c = getc(fd);
                while (c == '#')
                {
                    fscanf_s(fd, "%[^\n\r] ", b);
                    printf("%s\n", b);
                    c = getc(fd);
                }
                ungetc(c, fd);
            }
            fscanf_s(fd, "%d %d %d", &n, &m, &k);
            printf("%d rows  %d columns  max value= %d\n", m, n, k);
            c = getc(fd); // Skip the last whitespace

            numbytes = n * m * 3;
            image = (unsigned char *)malloc(numbytes);
            if (image == NULL)
            {
                printf("Memory allocation failed!\n");
                return NULL;
            }
            // Read and re-order as necessary
            for (i = m - 1; i >= 0; i--) for (j = 0; j < n; j++) // Important bug fix here!
            {
                image[(i*n + j) * 3 + 0] = getc(fd);
                image[(i*n + j) * 3 + 1] = getc(fd);
                image[(i*n + j) * 3 + 2] = getc(fd);
            }
        }
        else
        {
            printf("%s is not a PPM file!\n", filename);
            return NULL;
        }

    printf("read image\n");

    *height = m;
    *width = n;
    return image;
}