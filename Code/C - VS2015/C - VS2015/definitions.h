#pragma once

#define M_2_PI 6.28318530718
#define N 33554432 // 73 times -> ~180

typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;

struct complex
{
	float r;
	float i;
};

const int tab32[32] = { 0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30, 8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31 };