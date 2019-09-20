#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <omp.h>

typedef unsigned char pixel_t[3]; // colors [R, G ,B]

// picture resolution
static const int ImageWidth = 2000;
static const int ImageHeight = 2000;
static const double CxMin = -2.5;
static const double CxMax = 1.5;
static const double CyMin = -2.0;
static const double CyMax = 2.0;

const int IterationMax = 150;
const int MaxColorComponentValue = 1 << 8; //256
// const double Bailout = 2;                       // bail-out value
const double Circle_Radius = 4; //Bailout * Bailout; // circle radius

int main()
{

    double PixelWidth = (CxMax - CxMin) / ImageWidth;   /* scaled x coordinate of pixel (must be scaled to lie somewhere in the Mandelbrot X scale (-2.5, 1.5) */
    double PixelHeight = (CyMax - CyMin) / ImageHeight; /* scaled y coordinate of pixel (must be scaled to lie somewhere in the Mandelbrot Y scale (-2.0, 2.0) */

    //create array the size of the number of pixels in the image
    pixel_t *pixels = malloc(sizeof(pixel_t) * ImageHeight * ImageWidth);

    double run_time = omp_get_wtime();
    double Cy;
#pragma omp parallel for private(Cy) shared(pixels)
    for (int iY = 0; iY < ImageHeight; iY++)
    {
        Cy = CyMin + iY * PixelHeight;
        if (fabs(Cy) < PixelHeight / 2)
        {
            Cy = 0.0; // Main antenna
        }
        // #pragma omp parallel for
        for (int iX = 0; iX < ImageWidth; iX++)
        {
            double Cx = CxMin + iX * PixelWidth;
            double Zx = 0.0;
            double Zy = 0.0;
            double Zx2 = Zx * Zx;
            double Zy2 = Zy * Zy;

            int Iteration;
            for (Iteration = 0; Iteration < IterationMax && ((Zx2 + Zy2) < Circle_Radius); Iteration++)
            { //
                Zy = 2 * Zx * Zy + Cy;
                Zx = Zx2 - Zy2 + Cx;
                Zx2 = Zx * Zx;
                Zy2 = Zy * Zy;
            };

            if (Iteration == IterationMax)
            {
                //  interior of Mandelbrot set = black
                pixels[iY * ImageWidth + iX][0] = 0;
                pixels[iY * ImageWidth + iX][1] = 0;
                pixels[iY * ImageWidth + iX][2] = 0;
            }
            //
            else
            {
                pixels[iY * ImageWidth + iX][0] = ((double)(Iteration - log2(log2(sqrt(Zx2 + Zy2)))) / IterationMax) * MaxColorComponentValue;
                pixels[iY * ImageWidth + iX][1] = 0;
                pixels[iY * ImageWidth + iX][2] = 0;
            }
        }
    }
    run_time = omp_get_wtime() - run_time;

    printf("The time for Mandelbrot set is: %f \n", run_time);

    FILE *fp;
    fp = fopen("MandelbrotSet.ppm", "wb");
    fprintf(fp, "P6\n %s\n %d\n %d\n %d\n", "# no comment", ImageWidth, ImageHeight, MaxColorComponentValue);
    for (int iY = 0; iY < ImageHeight; iY++)
        for (int iX = 0; iX < ImageWidth; iX++)
            fwrite(pixels[iY * ImageWidth + iX], 1, sizeof(pixel_t), fp);
    fclose(fp);
    free(pixels);
    //  stop_timer ( );
    //
    //  printf("Elapsed time: %lf\n",elapsed_time ( ));
    return 0;
}