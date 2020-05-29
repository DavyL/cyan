#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include <cyan/color/color.h>
#include <cyan/image/image.h>
#include <cyan/image/load_png.h>
#include <cyan/image/complex.h>
#include <cyan/image/transforms.h>
#include <cyan/image/sampling.h>


//TODO : use this sampling function in FFT * to obtain odd and even elem 
int sub_sampling_step( void ** dst, void *  src, int number_elem, int step, size_t size, void (*assign)(void *, void *) ){
	
	if(dst == NULL){
		fprintf(stderr, "ERROR : sub_sampling_step: dst is a NULL pointer.\n");
		return -1;
	}

	if(*dst == NULL){
		fprintf(stdout, "sub_sampling_step : allocating memory ( allocated : %d, \ttotal size : %d ).\n",(number_elem/step)*size , number_elem * size );
		*dst = malloc( (number_elem/step)  * size );	
	}	
	int i = 0;
	for( i =0; i*step < number_elem; i++){
		fprintf(stdout, "sub_sampling_step : assigning data ( %d ). \n", i * step * size );
		assign( *dst + i*size, src + i * step * size );		
	}

	fprintf(stdout, "sub_sampling_step : leaving function. \n");

	return 0;
}

int sub_sampling_step_polar( complex_polar_t ** dst, complex_polar_t * src, int number_elem, int step){
	return  sub_sampling_step( (void **) dst, src, number_elem, step, sizeof(complex_polar_t), (void (*)( void *, void *) ) assign_complex_polar_ptr );
}

int sub_sampling_step_cart( complex_cart_t ** dst, complex_cart_t * src, int number_elem, int step){
	return  sub_sampling_step( (void **) dst, src, number_elem, step, sizeof(complex_cart_t), (void (*)( void *, void *) ) assign_complex_cart_ptr );
}

double window_rectangular( double x, int N ){
	if ( x >= 0.0f && x <= (double) N - 1.0f )
		return 1.0f;
	return 0.0f;
}
double window_welch( double x, int N ){
	return 1.0f - pow( (( x -0.5f * ((double) N - 1.0f ) ) / (0.5f * ( (double) N + 1.0f )  ) ), 2.0f ) ;
}
double window_triangular( double x, int N ){
	return  (2.0f / ((double) N - 1.0f) ) * ( (double) N / 2.0f - fabs( x - (float) N /2.0f ) );
}
double window_hann( double x, int N){
	return 0.5f*(1.0f - cos( 8.0f * atan(1.0f) * x / (double) N - 1.0f ) );
}
double window_hamming( double x, int N){
	return 0.54f - 0.46 * cos( 8.0f * atan(1.0f) * x / (double) N - 1.0f );
}
int double_array_through_window( double ** dst, double * src, int number_elem, double (*window)(double, int) ){
	return arb_array_through_window( (void **) dst, src, number_elem, 
					sizeof(double ),
					window,
					assign_double,
					mult_double );

}
int assign_double( void * dst, void * src ){
	*((double *) dst) = *((double *) src);
	return 0;
}
int mult_double( void * dst, double src ){
	*((double *) dst) *= src;
	return 0;
}
int cart_array_through_window( complex_cart_t ** dst, complex_cart_t * src, int number_elem, double (*window)(double, int) ){
	if(*dst == (complex_cart_t * ) NULL){
		*dst = malloc( number_elem * sizeof(complex_cart_t) );
		if(*dst == (complex_cart_t * )  NULL){
			fprintf(stderr, "array_through_window : Error allocating memory.\n");
			return -1;
		}
	}

	int i = 0;
	for(i = 0; i < number_elem; i++){
		(*dst)[i].real = src[i].real * window(i, number_elem);
		(*dst)[i].im = src[i].im * window(i, number_elem);
	}
}
int arb_array_through_window( void ** ret, void * array, int number_elem, size_t elem_size, double (*window)(double, int), int (*assign)(void *, void *), int (*mult)( void *, double) ){
	if(*ret == NULL){
		*ret = malloc( number_elem * elem_size );
		if(*ret == NULL){
			fprintf(stderr, "array_through_window : Error allocating memory.\n");
			return -1;
		}
	}

	int i = 0;
	for(i = 0; i < number_elem; i++){
		assign( *ret + i * elem_size, array + i * elem_size);
		mult( *ret + i *elem_size, window((double) i, number_elem) );
	}
	return 0;
}

int cart_array_through_window_arb( complex_cart_t ** ret, complex_cart_t * array, int N, double (*window)(double, int) ){
	return arb_array_through_window( (void **) ret, array, N, 
					sizeof(complex_cart_t ),
					window,
					assign_complex_cart_ptr_void,
					mult_scalar_complex_cart );

}
int double_array_through_real_fct( double ** ret, int N, double (*fct)(double, int) ){
	double  * ones =  malloc( N *sizeof(double) );
	int i = 0;
	for(i = 0; i < N; i++)
		ones[i] = 1.0f;
	int temp = arb_array_through_window( (void **) ret, ones, N, sizeof(double), fct, assign_double, mult_double); 
	free(ones);
	return temp;
}
//Return index number of max value in array
int max_index_in_double_array(double * array, int len){
	if(array == NULL){
		fprintf(stderr, "max_in_double_array : array is NULL\n");
		return -1;
	}

	int index_max = 0;
	double max = array[index_max]; //max is used to avoid retrieving data in array in each loop
	int i = 1;
	for( ; i < len; ++i){
		if(array[i] > max){
			max = array[i];
			index_max = i;
		}
	}

	return index_max;
}
int min_index_in_double_array(double * array, int len){
	if(array == NULL){
		fprintf(stderr, "min_in_double_array : array is NULL\n");
		return -1;
	}
	int index_min = 0;
	int min = array[0];
	int i = 1;
	for( ; i < len; ++i){
		if(array[i] < min){
			min = array[i];
			index_min = i;
		}
	}
	return index_min;
}

double max_val_in_double_array(double * array, int len){
	int index = max_index_in_double_array(array, len);
	if(index == -1){
		fprintf(stderr, "max_val_in_double_array : couldn't get max value, array is NULL\n");
		return 0.0f;
	}
	return array[index];
}
double min_val_in_double_array(double * array, int len){
	int index = min_index_in_double_array(array, len);
	if(index == -1){
		fprintf(stderr, "min_val_in_double_array : couldn't get min value, array is NULL\n");
		return 0.0f;
	}
	return array[index];
}

int normalize_and_scale_double_array(double * array, int len, double scale){
	if(array == NULL){
		fprintf(stderr, "normalize_and_scale_double_array : array is NULL\n");
		return -1;
	}
	int i = 0;
	double min = min_val_in_double_array(array, len);
	double normalizer = fabs( max_val_in_double_array(array, len) - min);
	
	if(normalizer == 0.0f){
		fprintf(stdout, "normalize_and_scale_double_array : array contains constant values, setting them to scale.\n");
		for(i = 0; i < len; i++){
			array[i] = scale;
		}
		return len;
	}
	add_double_array( array, len, - min);  //Maps all values to [0, max - min ]
	normalizer = scale / normalizer;

	return mult_double_array(array, len, normalizer);	//Maps all values to [0, scale]
}
//pw_mult_polar_array: pointwise multiplication of two polar arrays
int pw_mult_polar_array( complex_polar_t ** dst, complex_polar_t * x, complex_polar_t * y, int len){
	if(dst == NULL || x == NULL || y == NULL){
		fprintf(stdout, "pw_mult_polar_array: A NULL was given as an argument.\n");
		return -1;
	}
	if(*dst == NULL){
		*dst = malloc(len * sizeof(complex_polar_t));
		if(*dst == NULL){
			fprintf(stderr, "pw_mult_polar_array: couldn't allocate memory, len : %d.\n", len);
			return -1;
		}
	}
	int i = 0;
	for(i = 0; i < len; i++){
		(*dst)[i].power = x[i].power * y[i].power;
		(*dst)[i].phase = x[i].phase + y[i].phase;
	}
	return len;
}
int double_array_to_polar_power_array(complex_polar_t ** dst, double * src, int len){
	if(dst == NULL || src == NULL){
		fprintf(stdout, "double_array_to_polar_array : A NULL was given as an arguments.\n");
		return -1;
	}
	if(*dst == NULL){
		*dst = malloc(len * sizeof(complex_polar_t));
		if(*dst == NULL){
			fprintf(stderr, "double_array_to_polar_power_array : couldn't allocate memory : len : %d.\n", len);
			return -1;
		}
	}
	int i = 0;
	for( i = 0; i < len; i++){
		(*dst)[i].power = src[i];
		(*dst)[i].phase = 0.0f;
	}
	return 0;
}	
//mult_double_array: multiplication of an array by a fixed double
int mult_double_array(double * array, int len, double mult){
	if( array == NULL){
		fprintf(stdout, "normalize_double_array : array is NULL. \n");
		return -1;
	}
	int i = 0;
	for(i = 0; i < len; i++){
		array[i] *= mult;
	}
	return len;
}

int add_double_array(double * array, int len, double add){
	if( array == NULL){
		fprintf(stdout, "add_double_array : array is NULL\n");
		return -1;
	}
	int i = 0;
	for(i = 0; i < len; i++){
		array[i] += add;
	}
	return len;
}

int double_array_to_int_array(int ** dst, double * src, int len){
	if(dst == NULL || src == NULL){
		fprintf(stdout, "double_array_to_int_array : dst or src is NULL\n");
		return -1;
	}
	*dst = malloc(len * sizeof(int));
	if(*dst == NULL){
		fprintf(stdout, "double_array_to_int_array : allocation error\n");
		return -1;
	}
	
	int i = 0;
	for( ; i < len; i++){
		(*dst)[i] = (int) src[i];
	}
	return len;

}

int array_2d_through_window_arb( void *** dst, void ** src, int N, int M, size_t elem_size, 
	       			double (*window)(double, int), int (*mult)(void *, double)){
 	
	int i, j;
	*dst = (void **)  malloc( N * sizeof( void *) );

	for(i = 0; i < N; i++){
		(*dst)[i] = malloc( N * elem_size );
		if( (*dst)[i] == NULL){
			fprintf(stderr, "array_2d_through_window_arb : memory allocation failed.\n");
			return -1;
		}
		
		cart_array_through_window_arb( (*dst)[i], *(src + i), M, window);	
	}
	for(j = 0; j < M; j++){
		for( i = 0; i < N; i++){
			mult( (*dst)[j] + i * elem_size, window( (double) j, N) );
		}
	}

	return 0;
}

