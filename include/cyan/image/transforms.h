#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include <cyan/image/image.h>
#include <cyan/image/complex.h>

void * rotate_buffer( void * buffer, int N, size_t );
void * reverse_buffer( void * buffer, int N, size_t );
int transpose( void *** dst, void ** src, int N, size_t size );
int rotate_180( void *** dst, void ** src, int N, size_t size);

image_t * FT_image_Y(image_t * image, complex_cart_t ** (*transform)(complex_cart_t **, int, int) , image_t * (*convert_cart_to_image)( complex_cart_t **, int, int));
complex_cart_t ** image_to_cart(image_t * image, complex_cart_t (*xyz_to_cart)( double, double, double) );
image_t * cart_to_image(image_t * image, complex_cart_t ** array, int n, int m, double (*X_func)(complex_cart_t), double (*Y_func)(complex_cart_t), double (*Z_func)(complex_cart_t));
image_t * cart_to_power_Y_image(image_t * image, complex_cart_t ** array, int N, int M);

complex_cart_t ** lowpass_harsh(complex_cart_t ** array_cart, int N, int M);
complex_cart_t ** highpass_harsh(complex_cart_t ** array_cart, int N, int M);
complex_cart_t ** middle_band_pass( complex_cart_t ** array_cart, int N, int M);

complex_cart_t Y_to_cart(double X, double Y, double Z);
complex_cart_t X_Y_to_cart(double X, double Y, double Z);
complex_cart_t X_Y_Z_to_cart(double X, double Y, double Z);

image_t * cart_to_Y(complex_cart_t ** array, int N, int M);
image_t * cart_power_to_Y(complex_cart_t ** array, int N, int M);
image_t * cart_log_power_to_Y(complex_cart_t ** array, int N, int M);
image_t * cart_phase_to_Y(complex_cart_t ** array, int N, int M);

void * unitary_ft_polar(complex_polar_t * polar_array, int N);
void * unitary_ft_cart(complex_cart_t * cart_array, int N);

image_t * color2grey(image_t *);

float norm_2_sq_4_f( float, float , float , float );
float dot_product_4_f(	float, float, float, float, float, float, float, float );

int double_threshold( double * array, double T, int len);
int l2_distance(double * dst, double * x, double * y, int len);

int convol_loc( double * a, double * array, double * h, int i, int filter_len, int len );


#endif
