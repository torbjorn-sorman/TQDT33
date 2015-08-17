#include "mrfft.h"

/* Naive Fast Fourier Transform */
void mrfft(int dir, my_complex *x, my_complex *X, uint32_t N)
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

void mrfft2d(int dir, void(*fn)(int, my_complex*, my_complex*, uint32_t), my_complex **seq2d, uint32_t N)
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

        fn(dir, seq, out, N);

        for (col = 0; col < N; ++col)
            seq2d[row][col] = out[col];
    }

    for (col = 0; col < N; ++col)
    {
        for (row = 0; row < N; ++row)
            seq[row] = seq2d[row][col];

        fn(dir, seq, out, N);

        for (row = 0; row < N; ++row)
            seq2d[row][col] = out[row];
    }
    free(seq);
}

/* Naive Discrete Fourier Transform, essentially as per definition */
void naive_dft(my_complex *x, my_complex *X, uint32_t N)
{
    float real, img;
    my_complex y = { 0.0, 0.0 };
    float re, im;
    my_complex tmp = { 0.0, 0.0 };
    float theta = 1.0;
    float c1 = -M_2_PI / N;
    float c2 = 1.0;
    for (uint32_t k = 0; k < N; ++k)
    {
        real = 0.0;
        img = 0.0;
        c2 = c1 * k;
        for (uint32_t n = 0; n < N; ++n)
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