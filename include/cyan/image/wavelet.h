#ifndef WAVELET_H
#define WAVELET_H

#include <cyan/image/complex.h>

int D4_wavelet_forward( double ** dst_scale, double ** dst_wavelet, double * src, int len, int offset);
int D4_wavelet_backward( double ** dst, double * src_scale, double * src_wavelet, int len, int offset);

int D4_wavelet_forward_columns( double ** dst_scale, double ** dst_wavelet, double * src, int len_row, int len_col);
int D4_wavelet_forward_rows( double ** dst_scale, double ** dst_wavelet, double * src, int len_row, int len_col);

int D4_wavelet_2D( double ** sca_sca, double ** sca_wav, double ** wav_sca, double ** wav_wav, double * array, int len_row, int len_col);
int D4_wavelet_2D_CR( double ** sca_sca, double ** sca_wav, double ** wav_sca, double ** wav_wav, double * array, int len_row, int len_col);

int D4_wavelet_backward_rows( double ** dst, double * src_scale, double * src_wavelet, int len_row, int len_col);
int D4_wavelet_backward_columns( double ** dst, double * src_scale, double * src_wavelet, int len_row, int len_col);

int D4_wavelet_2D_backward( double ** dst, double * sca_sca, double * sca_wav, double * wav_sca, double * wav_wav, int len_row, int len_col);
int D4_wavelet_2D_backward_CR( double ** dst, double * sca_sca, double * sca_wav, double * wav_sca, double * wav_wav, int len_row, int len_col);

int compute_dual_filter(double ** filter_dst, double * filter_src, int filter_length);
int D4_wavelet_2D_convol( double ** sca_sca, double ** sca_wav, double ** wav_sca, double ** wav_wav, double * array, int len_row, int len_col);

int wavelet_2D_sampling(double ** dst, double * src, int src_rows, int src_cols);

#endif
