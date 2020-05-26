#ifndef FOURIER_H
#define FOURIER_H

#include <cyan/image/complex.h>
#include <cyan/image/transforms.h>

complex_polar_t * FFT_1D( complex_polar_t *, complex_polar_t *, int);
complex_polar_t * FFT_1D_reverse( complex_polar_t * fourier_polar, int n);

complex_cart_t * FFT_1D_cart_to_cart( complex_cart_t * f, int n);
complex_cart_t * FFT_1D_reverse_cart_to_cart( complex_cart_t * fft, int n);

complex_cart_t ** FFT_2D(complex_cart_t ** array_cart, int n, int m);
complex_cart_t ** FFT_2D_reverse(complex_cart_t ***, complex_cart_t ** ft_array_cart, int n, int m);




#endif

