#include <stdio.h>
#include <Windows.h>
#include <cmath>

#include "definitions.h"

int log2_32(uint32_t value)
{
	value |= value >> 1; value |= value >> 2; value |= value >> 4; value |= value >> 8; value |= value >> 16;
	return tab32[(uint32_t)(value * 0x07C4ACDD) >> 27];
}

uint32_t reverse(uint32_t x, uint32_t l)
{
	x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
	x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
	x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
	x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
	return((x >> 16) | (x << 16)) >> (32 - l);
}

void naive_dft(complex *x, complex *X);
void fft(complex *x, complex *X);
void printTime(LARGE_INTEGER tStart, LARGE_INTEGER tStop, LARGE_INTEGER freq);
void printResult(complex *c, int n, char *str, int verified);
int verify_impulse(complex *c, int size);
void compareComplex(complex *c1, complex *c2, float t1, float t2);
void test_p();

int main()
{
	LARGE_INTEGER freq, tStart, tStop;
	float time_FFT, time_NDFT;

	/* Get ticks per second */
	QueryPerformanceFrequency(&freq);

	/* Prep data */
	complex *impulse = (complex *)malloc(sizeof(complex)*N);
	complex *res_FFT = (complex *)malloc(sizeof(complex)*N);
	for (int i = 0; i < N; ++i)
	{
		impulse[i] = { 0.0, 0.0 };
		res_FFT[i] = { 0.0, 0.0 };
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
	time_FFT = (float)(tStop.QuadPart - tStart.QuadPart) * 1000.0 / (float)freq.QuadPart;
	/*
	printResult(res_FFT, N, "impulse", verify_impulse(res_FFT, N));
	compareComplex(res_FFT, res_NDFT, time_FFT, time_NDFT);
	*/
	free(impulse);
	free(res_FFT);
	getchar();
	return 0;
}

/* Naive Discrete Fourier Transform, essentially as per definition */
void naive_dft(complex *x, complex *X)
{
	float real, img;
	complex y = { 0.0, 0.0 };
	float re, im;
	complex tmp = { 0.0, 0.0 };
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

/* Fast Fourier Transform */
void fft(complex *x, complex *X)
{
	complex *tmp = (complex *)malloc(sizeof(complex)*N);
	complex *W = (complex *)malloc(sizeof(complex)*N);
	const uint32_t depth = log2_32(N);
	const uint32_t n_half = N / 2;
	float w_angle = -M_2_PI / N;
	float theta;
	uint32_t bit = 0;
	float u_re, u_im, l_re, l_im;
	uint32_t u, l, p;
	uint32_t dist, dist_2, offset;

	for (uint32_t n = 0; n < N; ++n)
	{
		tmp[n] = x[n];
		p = n;
		p = (((p & 0xaaaaaaaa) >> 1) | ((p & 0x55555555) << 1));
		p = (((p & 0xcccccccc) >> 2) | ((p & 0x33333333) << 2));
		p = (((p & 0xf0f0f0f0) >> 4) | ((p & 0x0f0f0f0f) << 4));
		p = (((p & 0xff00ff00) >> 8) | ((p & 0x00ff00ff) << 8));
		theta = w_angle * (((p >> 16) | (p << 16)) >> (32 - depth));
		W[n].r = cos(theta);
		W[n].i = sin(theta);
	}

	dist = N;
	for (uint32_t k = 0; k < depth; ++k)
	{
		bit = depth - 1 - k;
		dist_2 = dist;
		dist = dist >> 1;
		offset = 0;
		for (uint32_t n = 0; n < n_half; ++n)
		{
			offset += (n >= (dist + offset)) * dist_2;
			l = (n & ~(1 << bit)) + offset;
			u = l + dist;
			// Lower			
			p = (l >> bit);
			l_re = tmp[l].r + W[p].r * tmp[u].r + W[p].i * tmp[u].i;
			l_im = tmp[l].i + W[p].r * tmp[u].i + W[p].i * tmp[u].r;
			// Upper
			p = (u >> bit);
			u_re = tmp[l].r + W[p].r * tmp[u].r + W[p].i * tmp[u].i;
			u_im = tmp[l].i + W[p].r * tmp[u].i + W[p].i * tmp[u].r;
			// Insert
			tmp[l].r = l_re;
			tmp[l].i = l_im;
			tmp[u].r = u_re;
			tmp[u].i = u_im;
		}
	}
	for (int n = 0; n < N; ++n)
	{
		X[n] = tmp[reverse(n, depth)];
	}
	free(tmp);
	free(W);
}

void printTime(LARGE_INTEGER tStart, LARGE_INTEGER tStop, LARGE_INTEGER freq)
{
	printf("Time (ms): %f\n", (float)(tStop.QuadPart - tStart.QuadPart) * 1000.0 / (float)freq.QuadPart);
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

void compareComplex(complex *c1, complex *c2, float t1, float t2)
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
	printf("\nImproved %d times.", (int)(t2 / t1));
	printf("\nTheoretical limit %d times.", (int)(((float)N * N) / ((float)N * log2_32(N))));
}