#ifndef FIR_FILTER_H_
#define FIR_FILTER_H_
#include <stdint.h>

extern void highpassfilter_coefficient(float *filter, float f, float fs, int N);
extern void lowpassfilter_coefficient(float *filter, float f, float fs, int N);
extern void bandpassfilter_coefficient(float *filter, float fl, float fh, float fs, int N);
extern void convolution(float *data, float *filter, uint32_t num, int N);
extern void reverce_array(float *arr, uint32_t n_data);
#endif
