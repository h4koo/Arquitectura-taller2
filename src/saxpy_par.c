#include <stdio.h>
#include <omp.h>

#define ARRAY_SIZE 100000

void saxpy(int n, float a, float *restrict x, float *restrict y)
{
#pragma omp parallel for num_threads(6)
    for (int i = 0; i < n; ++i)
        y[i] = a * x[i] + y[i];
}

void fillArray(int size, float *array, float value)
{
#pragma omp parallel for num_threads(6)
    for (int i = 0; i < size; ++i)
    {
        array[i] = value;
    }
}

int main()
{

    float x[ARRAY_SIZE];
    float y[ARRAY_SIZE];

    double run_time;

    run_time = omp_get_wtime();

    fillArray(ARRAY_SIZE, x, 3);
    fillArray(ARRAY_SIZE, y, 2);

    saxpy(ARRAY_SIZE, 4.5, x, y);

    run_time = omp_get_wtime() - run_time;

    printf("The time for SAXPY parallel is: %f \n", run_time);

    return 0;
}