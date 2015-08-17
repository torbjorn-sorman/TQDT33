#include <Windows.h>
#include <limits>

#include "mrimage.h"

void imgToComplex(unsigned char *img, my_complex **com, uint32_t N)
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

void minRange(double *m, double *r, my_complex **com, uint32_t N)
{
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
/*
void realToImg(my_complex **com, unsigned char *img, uint32_t N)
{
    uint32_t px;
    double vmin, range;
    minRange(&vmin, &range, com, N);
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
*/

void fftToImg(my_complex **com, unsigned char *img, uint32_t N)
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

void fftShift(unsigned char *in, unsigned char *out, uint32_t N)
{
    uint32_t n2 = (N / 2) * 3;
    uint32_t n4 = (N / 4) * 3;
    uint32_t _3n4 = n2 + n4;
    for (uint32_t i = 0; i < n4; ++i)
    {
        out[i + _3n4] = in[i];
    }
    for (uint32_t i = n4; i < n2; ++i)
    {
        out[i + n4] = in[i];
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