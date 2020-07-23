#ifndef WAVELET_H
#define WAVELET_H

#include <cyan/image/complex.h>

int D4_wavelet_forward( double ** dst_scale, double ** dst_wavelet, double * src, int len, int offset);
int D4_wavelet_backward( double ** dst, double * src_scale, double * src_wavelet, int len, int offset);

int D4_wavelet_forward_columns( double ** dst_scale, double ** dst_wavelet, double * src, int len_row, int len_col);
int D4_wavelet_forward_rows( double ** dst_scale, double ** dst_wavelet, double * src, int len_row, int len_col);

int D4_wavelet_2D( double ** sca_sca, double ** sca_wav, double ** wav_sca, double ** wav_wav, double * array, int len_row, int len_col);
int D4_wavelet_2D_CR( double ** sca_sca, double ** sca_wav, double ** wav_sca, double ** wav_wav, double * array, int len_row, int len_col);
#endif
