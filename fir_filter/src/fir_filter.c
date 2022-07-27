#include "fir_filter.h"
#include <math.h>

#define PI 3.14159265358979323846

void highpassfilter_coefficient(float *filter, float f, float fs, int N)
{
    float   fc = 0;
    float   omega = 0;
    int middle = 0;

    middle = (int)(N / 2);
    fc = f / fs;
    omega = 2 * PI * fc;

    for (int i = 0; i < N; i++) {
        if (i == middle) {
            filter[i] = 1 - 2 * fc;
        }
        else {
            filter[i] = (0 - sinf(omega * (i - middle))) / (PI * (i - middle));
        }
    }
}

void lowpassfilter_coefficient(float *filter, float f, float fs, int N)
{
    float   fc = 0;
    float   omega = 0;
    int middle = 0;

    middle = (int)(N / 2);
    fc = f / fs;
    omega = 2 * PI * fc;

    for (int i = 0; i < N; i++) {
        if (i == middle) {
            filter[i] = 2 * fc;
        }
        else {
            filter[i] = sinf(omega * (i - middle)) / (PI * (i - middle));
        }
    }
}

void bandpassfilter_coefficient(float *filter, float fl, float fh, float fs, int N)
{
    float   flc = 0, omegal = 0;
    float   fhc = 0, omegah = 0;
    int middle = 0;

    middle = (int)(N / 2);
    flc = fl / fs;
    omegal = 2 * PI * flc;
    fhc = fh / fs;
    omegah = 2 * PI * fhc;

    for (int i = 0; i < N; i++) {
        if (i == middle) {
            filter[i] = 2 * fhc - 2 * flc;
        }
        else {
            filter[i] = \
                        sinf(omegah * (i - middle)) / (PI * (i - middle)) - \
                        sinf(omegal * (i - middle)) / (PI * (i - middle));
        }
    }
}

void convolution(float *data, float *filter, uint32_t num, int N)
{
    float   sum;

    for (uint32_t j = 0; j < num; j++) {
        sum = 0;
        for (uint32_t i = 0; i < N; i++) {
            sum += data[i + j] * filter[i];
        }
        data[j] = sum;
    }
}

void reverce_array(float *arr, uint32_t n_data)
{
    float   temp = 0.0f, sum = 0.0f, avg;

    for (uint32_t i = 0, j = n_data - 1; i < n_data / 2; i++, j--) {
        temp = arr[j];
        arr[j] = arr[i];
        arr[i] = temp;
        sum += arr[j] +arr[i];
    }

    avg = sum / n_data;

    for (uint32_t i = 0; i < n_data; i++) {
        arr[i] -= avg;
    }
}
