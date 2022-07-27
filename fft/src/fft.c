#include <math.h>
#include <stdint.h>
#include "fft.h"

static void swap(float *a, float *b);

#define _PI     3.14159265359f
#define _2PI    6.28318530718f
#define HANNING(IN, IDX, LEN) do { \
    (IN) = (IN) * ((1.0 - cos(2.0 * _PI * (float)(IDX) / (float)((LEN) - 1)))); \
} while (0)

void hanning_window(float *buff, uint32_t len)
{
    for (uint32_t i = 0; i < len; i++) {
        HANNING(buff[i], i, len);
    }
}

void gen_fft_table(float wn_FFT[], short br_FFT[], int N_FFT)
{
    int n_half;
    float   arg;

    arg = _2PI / N_FFT;
    for (int i = 0; i < ((N_FFT * 3) >> 2); i++) {
        wn_FFT[i] = cosf(arg * i);
    }

    n_half = N_FFT >> 1;
    br_FFT[0] = 0;
    for (int ne = 1; ne < N_FFT; ne = ne << 1) {
        for (int jp = 0; jp < ne; jp++) {
            br_FFT[jp + ne] = br_FFT[jp] + n_half;
        }
        n_half = n_half >> 1;
    }
}

void fft(float xR[], float xI[], float wn_FFT[], short br_FFT[], int N_FFT)
{
    float   xtmpR, xtmpI;
    int jnh, jxC, jxS, n_half, n_half2;

    n_half = N_FFT >> 1;

    for (int ne = 1; ne < N_FFT; ne = ne << 1) {
        n_half2 = n_half << 1;
        for (int k = 0; k < N_FFT; k = k + n_half2) {
            jxC = 0;
            jxS = N_FFT >> 2;
            for (int j = k; j < (k + n_half); j++) {
                jnh = j + n_half;

                xtmpR = xR[j];
                xtmpI = xI[j];
                xR[j] = xtmpR + xR[jnh];
                xI[j] = xtmpI + xI[jnh];
                xtmpR = xtmpR - xR[jnh];
                xtmpI = xtmpI - xI[jnh];
                xR[jnh] = xtmpR * wn_FFT[jxC] - xtmpI * wn_FFT[jxS];
                xI[jnh] = xtmpR * wn_FFT[jxS] + xtmpI * wn_FFT[jxC];

                jxC = jxC + ne;
                jxS = jxS + ne;
            }
        }
        n_half = n_half >> 1;
    }

    for (int j = 0; j < N_FFT; j++)  {
        if (j < br_FFT[j]) {
            swap(&xR[j], &xR[br_FFT[j]]);
            swap(&xI[j], &xI[br_FFT[j]]);
        }
    }
}

static void swap(float *a, float *b)
{
    float tmp;

    tmp = *a;
    *a = *b;
    *b = tmp;
}
