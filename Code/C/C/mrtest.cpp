#include <Windows.h>
#include <limits>

#include "mrtest.h"
#include "mrfft.h"
#include "mrimage.h"

double run_test(void(*fn)(int, my_complex*, my_complex*, uint32_t), my_complex *in, my_complex *out, uint32_t N)
{
    LARGE_INTEGER freq, tStart, tStop;
    double m = DBL_MAX;
    QueryPerformanceFrequency(&freq);
    for (int i = 0; i < 10; ++i)
    {
        QueryPerformanceCounter(&tStart);
        fn(FORWARD_FFT, in, out, N);
        QueryPerformanceCounter(&tStop);
        m = min(m, (double)(tStop.QuadPart - tStart.QuadPart) * 1000.0 / (float)freq.QuadPart);
    }
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

double run_imgTest(char* filename, void(*fn)(int, my_complex*, my_complex*, uint32_t), uint32_t N)
{
    int n, m;
    /* Read image to memory */
    unsigned char *image = readppm(filename, &n, &m);
    if (n != N || m != N)
    {
        printf("Image size not square and pw of 2.\n");
        getchar();
        return -1.0;
    }
    /* Allocate resources */
    unsigned char *imImage = (unsigned char *)malloc(sizeof(unsigned char) * N * N * 3);
    unsigned char *imImage2 = (unsigned char *)malloc(sizeof(unsigned char) * N * N * 3);
    my_complex **cxImage = (my_complex **)malloc(sizeof(my_complex) * N);
    for (int i = 0; i < N; ++i)
    {
        cxImage[i] = (my_complex *)malloc(sizeof(my_complex) * N);
    }

    /* Set real values from image values. 
     * Store the real-value version of the image.
     */
    imgToComplex(image, cxImage, N);
    fftToImg(cxImage, imImage, N);
    printf("Write img00-grey.ppm\n");
    writeppm("img00-grey.ppm", N, N, imImage);    

    /* Run 2D FFT on complex values.
     * Map absolute values of complex to pixels and store to file.
     */
    mrfft2d(FORWARD_FFT, mrfft, cxImage, N);
    fftToImg(cxImage, imImage, N);
    //fftShift(imImage, imImage2, N);
    printf("Write img01-magnitude.ppm\n");
    writeppm("img01-magnitude.ppm", N, N, imImage);

    /* Run inverse 2D FFT on complex values */
    mrfft2d(INVERSE_FFT, mrfft, cxImage, N);
    fftToImg(cxImage, imImage, N);
    printf("Write img02-fftToImage.ppm\n");
    writeppm("img02-fftToImage.ppm", N, N, imImage);

    // Free all resources...
    free(image);
    for (int i = 0; i < N; ++i)
        free(cxImage[i]);    
    free(cxImage);
    free(imImage);
    free(imImage2);

    printf("Done...\n");
}

void fft_ref(int dir, my_complex *in, my_complex *out, uint32_t N)
{
    kiss_fft_cfg cfg = kiss_fft_alloc(N, (dir == -1), 0, 0);
    kiss_fft(cfg, in, out);
    free(cfg);
}