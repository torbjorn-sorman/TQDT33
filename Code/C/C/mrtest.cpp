#include <Windows.h>
#include <limits>

#include "mrtest.h"
#include "mrfft.h"
#include "mrimage.h"

#define DBG_PRINT

int run_simpleTest(uint32_t N)
{
    uint32_t i;
    my_complex *in = (my_complex *)malloc(sizeof(my_complex) * N);
    my_complex *mr_out = (my_complex *)malloc(sizeof(my_complex) * N);
    my_complex *re_out = (my_complex *)malloc(sizeof(my_complex) * N);
    double margin = 0.00001;
    // constant
    for (i = 0; i < N; ++i)
    {
        in[i].r = 0.5;//sin(M_2_PI * (((float)i) / N));
        in[i].i = 0.f;
        //printf("%f\n", in[i].r);
    }
    printf("\n");
    fft_ref(FORWARD_FFT, in, re_out, N);
    mrfft(FORWARD_FFT, in, mr_out, N);
    for (i = 0; i < N; ++i)
    {
        printf("{%f, %f} {%f, %f}\n", mr_out[i].r, mr_out[i].i, re_out[i].r, re_out[i].i);
    }
    for (i = 0; i < N; ++i)
    {
        if (abs(re_out[i].r - mr_out[i].r) > margin || abs(re_out[i].i - mr_out[i].i) > margin)
        {
            printf("i: %u\nValues are:\t%f, %f\nShould be:\t%f, %f\n", i, mr_out[i].r, mr_out[i].i, re_out[i].r, re_out[i].i);
            return 0;
        }
    }
    return 1;
}

double run_test(void(*fn)(int, my_complex*, my_complex*, uint32_t), uint32_t N)
{
    LARGE_INTEGER freq, tStart, tStop;
    my_complex *in = (my_complex *)malloc(sizeof(my_complex) * N);
    my_complex *out = (my_complex *)malloc(sizeof(my_complex) * N);
    double m = DBL_MAX;
    uint32_t i;
    for (i = 0; i < N; ++i)
    {
        in[i].r = 0.f;//sin(M_2_PI * ((float)i / N));
        in[i].i = 0.f;
    }
    in[1].r = 1.f;
    QueryPerformanceFrequency(&freq);
    for (int i = 0; i < 10; ++i)
    {
        QueryPerformanceCounter(&tStart);
        fn(FORWARD_FFT, in, out, N);
        QueryPerformanceCounter(&tStop);
        m = min(m, (double)(tStop.QuadPart - tStart.QuadPart) * 1000.0 / (float)freq.QuadPart);
    }
#ifdef DBG_PRINT
    for (i = 0; i < N; ++i)
    {
        printf("{%.2f,\t%.2f}\n", out[i].r, out[i].i);
    }
#endif
    return m;
}

double run_2dtest(void(*fn)(int, my_complex*, my_complex*, uint32_t), my_complex **in, uint32_t N)
{
    LARGE_INTEGER freq, tStart, tStop;
    double m = DBL_MAX;
    QueryPerformanceFrequency(&freq);
    for (int i = 0; i < 10; ++i)
    {
        QueryPerformanceCounter(&tStart);
        mrfft2d(FORWARD_FFT, fn, in, N);
        QueryPerformanceCounter(&tStop);
        m = min(m, (double)(tStop.QuadPart - tStart.QuadPart) * 1000.0 / (float)freq.QuadPart);
    }
    return (double)(tStop.QuadPart - tStart.QuadPart) * 1000.0 / (float)freq.QuadPart;
}

int run_fbtest(void(*fn)(int, my_complex*, my_complex*, uint32_t), uint32_t N)
{
    uint32_t i;
    my_complex *in = (my_complex *)malloc(sizeof(my_complex) * N);
    my_complex *out = (my_complex *)malloc(sizeof(my_complex) * N);
    float *reference = (float *)malloc(sizeof(float)*N);
    float margin = 0.00001f;
    for (i = 0; i < N; ++i)
    {
        reference[i] = 0;//sin(M_2_PI * ((float)i / N));
        in[i].r = reference[i];
        in[i].i = 0.f;
    }
    in[1].r = reference[1] = 1.f;
#ifdef DBG_PRINT
    printf("\nInitial:\n");
    for (i = 0; i < N; ++i)
        printf("{ %f,\t%fi }\n", in[i].r, in[i].i);
#endif
    fn(FORWARD_FFT, in, out, N);
#ifdef DBG_PRINT
    printf("\nFFT:\n");
    for (i = 0; i < N; ++i)
        printf("{ %f,\t%fi }\n", out[i].r, out[i].i);
#endif
    fn(INVERSE_FFT, out, in, N);
#ifdef DBG_PRINT
    printf("\nInverse FFT:\n");
    for (i = 0; i < N; ++i)
        printf("{ %f,\t%fi }\n", in[i].r, in[i].i);
#endif
    for (i = 0; i < N; ++i)
    {
        if (abs(reference[i] - in[i].r) > margin)
            return 0;
    }
    return 1;
}

int run_fft2dinvtest(void(*fn)(int, my_complex*, my_complex*, uint32_t), uint32_t N)
{
    uint32_t i, x, y;
    double margin;
    my_complex **seq2d, **ref;
    char *format = "{%.2f, %.2f} ";
    seq2d = (my_complex **)malloc(sizeof(my_complex) * N * N);
    ref = (my_complex **)malloc(sizeof(my_complex) * N * N);
    for (i = 0; i < N; ++i)
    {
        seq2d[i] = (my_complex *)malloc(sizeof(my_complex) * N);
        ref[i] = (my_complex *)malloc(sizeof(my_complex) * N);
    }
    for (y = 0; y < N; ++y)
    {
        for (x = 0; x < N; ++x)
        {
            seq2d[y][x].r = ref[y][x].r = 0.f;
            seq2d[y][x].i = ref[y][x].i = 0.f;
        }
    }
    seq2d[1][1].r = ref[1][1].r = 1.f;
#ifdef DBG_PRINT
    for (y = 0; y < N; ++y)
    {
        for (x = 0; x < N; ++x)
        {
            printf(format, seq2d[y][x].r, seq2d[y][x].i);
        }
        printf("\n");
    }
    printf("\n");
#endif
    mrfft2d(FORWARD_FFT, fn, seq2d, N);
#ifdef DBG_PRINT
    for (y = 0; y < N; ++y)
    {
        for (x = 0; x < N; ++x)
        {
            printf(format, seq2d[y][x].r, seq2d[y][x].i);
        }
        printf("\n");
    }
    printf("\n");
#endif
    mrfft2d(INVERSE_FFT, fn, seq2d, N);
#ifdef DBG_PRINT
    for (y = 0; y < N; ++y)
    {
        for (x = 0; x < N; ++x)
        {
            printf(format, seq2d[y][x].r, seq2d[y][x].i);
        }
        printf("\n");
    }
    printf("\n");
#endif
    margin = 0.00001;
    for (y = 0; y < N; ++y)
    {
        for (x = 0; x < N; ++x)
        {
            if (abs(seq2d[y][x].r - ref[y][x].r) > margin)
                return 0;
            if (abs(seq2d[y][x].i - ref[y][x].i) > margin)
                return 0;
        }
    }
    return 1;
}

int run_imgTest(void(*fn)(int, my_complex*, my_complex*, uint32_t), uint32_t N)
{
    int res, n, m;
    uint32_t i;
    char filename[13];
    unsigned char *image, *imImage, *imImage2, *greyImage;
    my_complex **cxImage;

    sprintf_s(filename, 13, "lena_%u.ppm", N);
    /* Read image to memory */
    image = readppm(filename, &n, &m);
    if (n != N || m != N)
    {
        printf("Image size not square and pw of 2.\n");
        getchar();
        return 0;
    }
    /* Allocate resources */
    greyImage = (unsigned char *)malloc(sizeof(unsigned char) * N * N * 3);
    imImage = (unsigned char *)malloc(sizeof(unsigned char) * N * N * 3);
    imImage2 = (unsigned char *)malloc(sizeof(unsigned char) * N * N * 3);
    cxImage = (my_complex **)malloc(sizeof(my_complex) * N);
    for (i = 0; i < N; ++i)
    {
        cxImage[i] = (my_complex *)malloc(sizeof(my_complex) * N);
    }

    /* Set real values from image values.
     * Store the real-value version of the image.
     */
    imgToComplex(image, cxImage, N);
    fftToImg(cxImage, greyImage, N);
    printf("Write img00-grey.ppm\n");
    writeppm("img00-grey.ppm", N, N, greyImage);

    /* Run 2D FFT on complex values.
     * Map absolute values of complex to pixels and store to file.
     */
    mrfft2d(FORWARD_FFT, mrfft, cxImage, N);
    //fftToImg(cxImage, imImage, N);
    fftToMagnitudeImg(cxImage, imImage, N);
    fftShift(imImage, imImage2, N);
    //fftShift(imImage, imImage2, N);
    printf("Write img01-magnitude.ppm\n");
    writeppm("img01-magnitude.ppm", N, N, imImage2);
    //writeppm("img02-magnitude.ppm", N, N, imImage2);

    /* Run inverse 2D FFT on complex values */
    mrfft2d(INVERSE_FFT, mrfft, cxImage, N);
    fftToImg(cxImage, imImage, N);
    printf("Write img02-fftToImage.ppm\n");
    writeppm("img02-fftToImage.ppm", N, N, imImage);

    res = 1;
    for (i = 0; i < N * N * 3; ++i)
    {
        if (abs(greyImage[i] - imImage[i]) > 1)
        {
            printf("\nAt: %u\nIs: %u\nShould be: %u\n", i, imImage[i], greyImage[i]);
            res = 0;
            break;
        }
    }
    // Free all resources...
    free(image);
    for (i = 0; i < N; ++i)
        free(cxImage[i]);
    free(cxImage);
    free(imImage);
    free(imImage2);
    return res;
}


int run_fft2dTest(void(*fn)(int, my_complex*, my_complex*, uint32_t), uint32_t N)
{
    int res, n, m;
    uint32_t i;
    char filename[13];
    unsigned char *image, *imImage;
    my_complex **cxImage;

    sprintf_s(filename, 13, "lena_%u.ppm", N);
    /* Read image to memory */
    image = readppm(filename, &n, &m);
    if (n != N || m != N)
    {
        printf("Image size not square and pw of 2.\n");
        getchar();
        return 0;
    }
    /* Allocate resources */
    imImage = (unsigned char *)malloc(sizeof(unsigned char) * N * N * 3);
    cxImage = (my_complex **)malloc(sizeof(my_complex) * N);
    for (i = 0; i < N; ++i)
    {
        cxImage[i] = (my_complex *)malloc(sizeof(my_complex) * N);
    }

    /* Set real values from image values.
    * Store the real-value version of the image.
    */
    imgToComplex(image, cxImage, N);
    fftToImg(cxImage, image, N);
    mrfft2d(FORWARD_FFT, mrfft, cxImage, N);
    mrfft2d(INVERSE_FFT, mrfft, cxImage, N);
    fftToImg(cxImage, imImage, N);
    writeppm("img-out.ppm", N, N, imImage);

    res = 1;
    for (i = 0; i < N * N * 3; ++i)
    {
        if (abs(image[i] - imImage[i]) > 1)
        {
            res = 0;
            break;
        }
    }
    // Free all resources...
    free(image);
    for (i = 0; i < N; ++i)
        free(cxImage[i]);
    free(cxImage);
    free(imImage);
    return res;
}

void fft_ref(int dir, my_complex *in, my_complex *out, uint32_t N)
{
    kiss_fft_cfg cfg = kiss_fft_alloc(N, (dir == INVERSE_FFT), NULL, NULL);
    kiss_fft(cfg, in, out);
    free(cfg);
}