#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cyan/image/image.h>
#include <cyan/image/fourier.h>
#include <cyan/image/sampling.h>
#include <cyan/image/transforms.h>

#include <cyan/color/color.h>

#include <cyan_fileio/save_ppm.h>

int D4_wavelet_forward( double ** dst_scale, double ** dst_wavelet, double * src, int len);
int D4_wavelet_backward( double ** dst, double * src_scale, double * src_wavelet, int len);
int double_threshold( double * array, double T, int len);
int l2_distance(double * dst, double * x, double * y, int len);
int convol_loc( double * a, double * array, double * h, int i, int filter_len, int len );

int main( int argc, char** argv, char* envv ) {

	int len = 512;
	int i = 0;
	int scale = 255;
	double array[len];

	for(i=0; i < len; i++){
		array[i] = cos( 20*2 * M_PI * (double) i / (double) len ); 
	}
	
	//Wavelet coefs
	double * scale_array = NULL;
	double * wavelet_array = NULL;
	D4_wavelet_forward(&scale_array, &wavelet_array, array, len);


	//Non-linear approx
	double T = 0.2f;
	double_threshold(scale_array, T, len/2);
	double_threshold(wavelet_array, T, len/2);

	double * reversed_array = NULL;
	D4_wavelet_backward(&reversed_array, scale_array, wavelet_array, len/2);

	free(scale_array);
	free(wavelet_array);

	double L2_diff = 0;
	l2_distance(&L2_diff, array, reversed_array, len);
	fprintf(stdout, "L2 difference of array and reversed array is : %f.\n", L2_diff);

	int * int_array = NULL;
	int * int_reversed_array = NULL;

	normalize_and_scale_double_array(array, len, (double) scale);
	double_array_to_int_array( &int_array, array, len);

	normalize_and_scale_double_array(reversed_array, len, (double) scale);	
	double_array_to_int_array( &int_reversed_array, reversed_array, len);
	free(reversed_array);

	image_t * image = image_new_empty(len, scale);
	
	color_t color;
	color_assign( &color, 1.0f, 1.0f, 1.0f, 0);

	image_set_color_to_coordinates(image, color, int_array, 1);
	free(int_array);

	color_assign( &color, 0.7f, 0.3f, 0.8f, 0);

	image_set_color_to_coordinates(image, color, int_reversed_array, 1);
	free(int_reversed_array);

	image_save_ppm(image, "image.ppm");
	
	image_free(image);

	return 0;
}

//Computes the Daubechies 4 forward wavelet transform with as in input, src, a double array of length len
//Output is stored in dst_*, which are arrays of length len / 2
int D4_wavelet_forward( double ** dst_scale, double ** dst_wavelet, double * src, int len){

	if(dst_scale == NULL || dst_wavelet == NULL || src == NULL){
		fprintf(stderr, "Error : D4_wavelet_forward() A NULL argument was passed.\n");
		return -1;
	}
	if(len % 2){
		fprintf(stderr, "Error : D4_wavelet_forward() len is an odd number.\n");
		return -1;
	}
	if(*dst_scale == NULL){
		*dst_scale = malloc((len/2) * sizeof(double));
	}
	if(*dst_wavelet == NULL){
		*dst_wavelet = malloc((len/2) * sizeof(double));
	}
	if(*dst_wavelet == NULL || *dst_scale == NULL){
		fprintf(stderr, "Error: D4_wavelet_forward() Error allocating dst pointers.\n");
		return -1;
	}

	//Loading Daubechies4 low pass filter (scaling function)
	const int filter_length = 4;
	double h[] = { 0.482962913f, 0.836516303f, 0.224143868f, -0.129409522f}; 
	
	//Loading the associated high pass filter (wavelet function)
	//TODO : implement general computation of g from h
	double g[] = { -0.12940952f, -0.22414387f, 0.8365163f, -0.48296291f};

	int i;	
	
	double temp = 0.0f;
	//Doing finite support convolutions w.r.t the associated filter
	for(i = 0; i < len / 2 ; i++){
		convol_loc( &temp, src, h, 2*i, filter_length, len ); 
		(*dst_scale)[i] = temp;
		convol_loc( &temp, src, g, 2*i, filter_length, len ); 
		(*dst_wavelet)[i] = temp;
	}

	return 0;
}

//Sets to 0 every value which is absolutely smaller than T
//Returns -1 if there is an error, otherwise it returns the number of values kept in the array
int double_threshold( double * array, double T, int len){
	if(array == NULL){
		fprintf(stderr, "Invalid argument : double_threshold() arrray is a NULL pointer.\n");
		return -1;
	}
	int i = 0;
	for( i = len - 1; i >= 0; i--){
		if(fabs(array[i]) < T){
			array[i] = 0.0f;
			len--;
		}
	}
	return len;
}
//Takes as an input scale and wavelet arrays given by forward Daubechies 4 wavelet transform of length len
//stores in dst the resulting backward transform of length 2 * len	
int D4_wavelet_backward( double ** dst, double * src_scale, double * src_wavelet, int len){
	if(dst == NULL || src_scale == NULL || src_wavelet == NULL){
		fprintf(stderr, "Invalid argument : D4_wavelet_backward() One of the arguments is a NULL pointer.\n");
		return -1;
	}
	if(len < 2){
		fprintf(stderr, "Invalid argument : D4_wavelet_backward() len is than 2.\n");
		return -1;
	}
	if(*dst == NULL){
		*dst = malloc(2 * len * sizeof(double));
		if(*dst == NULL){
			fprintf(stderr, "Error: D4_wavelet_backward : Couldn't allocate memory.\n");
			return -1;
		}
	}
	double conj_h[] = { 0.224143868f, 0.8365163f, 0.482962913f, -0.12940952f };
	double conj_g[] = {  -0.129409522f, -0.48296291f, 0.836516303f, -0.22414387f };

	int i = 0;
	int j = 0;
	(*dst)[j++] = src_scale[len - 1] * conj_h[0] + src_wavelet[len - 1] * conj_h[1] + src_scale[0]*conj_h[2] + src_wavelet[0]*conj_h[3];
	(*dst)[j++] = src_scale[len-1] * conj_g[0] + src_wavelet[len - 1] * conj_g[1] + src_scale[0]*conj_g[2] + src_wavelet[0]*conj_g[3];
	for(i = 0; i < len - 1; i++){
		(*dst)[j++] = src_scale[i] * conj_h[0] + src_wavelet[i] * conj_h[1] + src_scale[i+1]*conj_h[2] + src_wavelet[ i + 1]*conj_h[3];
		(*dst)[j++] = src_scale[i] * conj_g[0] + src_wavelet[i] * conj_g[1] + src_scale[i+1]*conj_g[2] + src_wavelet[ i + 1]*conj_g[3];
	}

	return 0;
}
//Compute the l2 distance of two arrays
int l2_distance(double * dst, double * x, double * y, int len){
	if(dst == NULL || x == NULL || y == NULL){
		fprintf(stderr, "Invalid argument: l2_distance() a NULL pointer was passed as an argument.\n");
		return -1;
	}
	
	*dst = 0.0f;
	int i = 0;	
	for(i = 0; i < len; i++)
		*dst += fabs(x[i]*x[i] - y[i]*y[i]);
	*dst = sqrt(*dst);
	return 0;


}
	
//Stores in a the convolution of array (of length n) 
//at step i with a filter h of length filter_len
int convol_loc( double * a, double * array, double * h, int i, int filter_len, int len ){
	if( a == NULL || array == NULL || h == NULL){
		fprintf(stderr, "ERR : wavelet1D : convol_loc() : a NULL pointer was passed as an argument.\n");
		return -1;
	}
	*a = 0.0f;

	int k = 0;
	int j = 0;
	if(filter_len % 2 != 0)
		filter_len++;

	for( k = ( 0 > i - filter_len / 2 ? 0 : i) ; k < ( len < i + filter_len ? len : i + filter_len); k++, j++){
		*a += array[k]*h[j];
	}
	
	return 0;


}

int circ_convol( double * dst, double * x, double * y, int len ){
	
	//Compute FFT of x and y
	complex_polar_t * x_polar = NULL;
	double_array_to_polar_power_array(&x_polar, x, len);	
	
	complex_polar_t * y_polar = NULL;
	double_array_to_polar_power_array(&y_polar, y, len);	
	
	complex_polar_t * ft_x_polar = NULL;
	ft_x_polar = FFT_1D( x_polar, NULL, len);
	unitary_ft_polar( ft_x_polar, len);
	
	complex_polar_t * ft_y_polar = NULL;
	ft_y_polar = FFT_1D( y_polar, NULL, len);
	unitary_ft_polar( ft_y_polar, len);
	
	//Compute ft_dst term by term product of ft_x and ft_y

	complex_polar_t * ft_dst = NULL;
	pw_mult_polar_array( &ft_dst, ft_x_polar, ft_y_polar, len);

	//Compute in dst the reverse fourier transform of ft_dst
	complex_polar_t * polar_dst = NULL;

	polar_dst = FFT_1D_reverse(ft_dst, 8);

	
	
	//Checking if polar_dst is real valued
	int i = 0;
	for(i=0; i < len; i++){
		fprintf(stdout,"%f\n", polar_dst[i].phase);
	}
//rotate_buffer( ft, len, sizeof(complex_polar_t));


	return 0;

}


