#ifndef SAMPLING_H
#define SAMPLING_H

#include "complex.h"

int sub_sampling_step( void **, void *, int, int, size_t, void (*)(void *, void * ) );
int sub_sampling_step_polar(complex_polar_t **, complex_polar_t *, int, int);
int sub_sampling_step_cart(complex_cart_t **, complex_cart_t *, int, int);

double window_rectangular( double x, int N );
double window_welch( double x, int N );
double window_triangular( double x, int N );
double window_hann( double x, int N);
double window_hamming( double x, int N);

int double_array_through_window( double ** dst, double * src, int number_elem, double (*window)(double, int) );
int assign_double( void * dst, void * src );
int mult_double( void * dst, double src );

int cart_array_through_window( complex_cart_t ** ret, complex_cart_t * array, int number_elem, double (*window)(double, int) );
int arb_array_through_window( void ** ret, void * array, int number_elem, size_t elem_size, double (*window)(double, int), int (*assign)(void *, void *), int (*mult)( void *, double) );
int cart_array_through_window_arb( complex_cart_t ** ret, complex_cart_t * array, int N, double (*)(double, int) );

int double_array_through_real_fct( double ** ret, int N, double (*fct)(double, int) );

int max_index_in_double_array(double * array, int len);
int min_index_in_double_array(double * array, int len);
double max_val_in_double_array(double * array, int len);
double min_val_in_double_array(double * array, int len);

int pw_mult_polar_array( complex_polar_t ** dst, complex_polar_t * x, complex_polar_t * y, int len);
int double_array_to_polar_power_array(complex_polar_t ** dst, double * src, int len);
int mult_double_array(double * array, int len, double mult);
int add_double_array(double * array, int len, double add);
int normalize_and_scale_double_array(double * array, int len, double scale);
int double_array_to_int_array(int ** dst, double * src, int len);
int double_array_to_polar_array( complex_polar_t ** dst, double * src, int len);

int array_2d_through_window_arb( void *** dst, void ** src, int N, int M, size_t elem_size, 
	       			double (*window)(double, int), int (*mult)(void *, double));


#endif

