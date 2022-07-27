#ifndef __FFT_H_
#define __FFT_H_
extern void gen_fft_table(float wn_FFT[], short br_FFT[], int N_FFT);
extern void fft(float xR[], float xI[], float wn_FFT[], short br_FFT[], int N_FFT);
extern void hanning_window(float *buff, uint32_t len);
#endif
