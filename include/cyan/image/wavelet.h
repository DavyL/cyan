#ifndef WAVELET_H
#define WAVELET_H

#include <cyan/image/complex.h>

int D4_wavelet_forward( double ** dst_scale, double ** dst_wavelet, double * src, int len);
int D4_wavelet_backward( double ** dst, double * src_scale, double * src_wavelet, int len);

#endif
