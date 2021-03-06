#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cyan/image/image.h>
#include <cyan/image/fourier.h>
#include <cyan/image/sampling.h>
#include <cyan/image/transforms.h>
#include <cyan/image/wavelet.h>

#include <cyan/color/color.h>

#include <cyan_fileio/save_ppm.h>

int main( int argc, char** argv, char* envv ) {

	int len = 512;
	int i = 0;
	int scale = 255;
	int offset = 1;
	double array[len];

	for(i=0; i < len; i++){
		array[i] = cos( 20*2 * M_PI * (double) i / (double) len ); 
	}
	
	//Wavelet coefs
	double * scale_array = NULL;
	double * wavelet_array = NULL;
	D4_wavelet_forward(&scale_array, &wavelet_array, array, len, offset);

	double norm_arr = 0;
	double norm_wav = 0;
	double norm_sca = 0;
	l2_distance(&norm_arr, array, NULL, len);
	l2_distance(&norm_wav, wavelet_array, NULL, len/2);
	l2_distance(&norm_sca, scale_array, NULL, len/2);
	fprintf(stdout, "array norm : %f\t, scale norm : %f\t, wavelet norm : %f\n", norm_arr, norm_sca, norm_wav); 



	//Non-linear approx
	double T = 0.2;
	double_threshold(scale_array, T, len/2);
	double_threshold(wavelet_array, T, len/2);

	double * reversed_array = NULL;
	D4_wavelet_backward(&reversed_array, scale_array, wavelet_array, len/2, offset);

	double L2_diff = 0;
	double norm_rev = 0;
	l2_distance(&norm_rev, reversed_array, NULL, len);
	fprintf(stdout, "reversed array norm : %f\n", norm_rev);
	l2_distance(&L2_diff, array, reversed_array, len);
	fprintf(stdout, "L2 difference of array and reversed array is : %f.\n", L2_diff);

	int * int_array = NULL;
	int * int_reversed_array = NULL;
	int * int_scale_array = NULL;
	int * int_wavelet_array = NULL;

	normalize_and_scale_double_array(array, len, (double) scale);
	double_array_to_int_array( &int_array, array, len);

	normalize_and_scale_double_array(reversed_array, len, (double) scale);	
	double_array_to_int_array( &int_reversed_array, reversed_array, len);

	normalize_and_scale_double_array( scale_array, len/2, (double) scale);
	double_array_to_int_array(&int_scale_array, scale_array, len/2);

	normalize_and_scale_double_array( wavelet_array, len/2, (double) scale);
	double_array_to_int_array(&int_wavelet_array, wavelet_array, len/2);

	image_t * image = image_new_empty(len, scale);
	image_t * image_wavelet = image_new_empty(len/2, scale);
	image_t * image_scale = image_new_empty(len/2, scale);
	
	color_t color;
	color_assign( &color, 1.0f, 1.0f, 1.0f, 0);

	image_set_color_to_coordinates(image, color, int_array, 1);
	free(int_array);

	color_assign( &color, 0.7f, 0.3f, 0.8f, 0);

	image_set_color_to_coordinates(image, color, int_reversed_array, 1);
	free(int_reversed_array);

	image_set_color_to_coordinates(image_wavelet, color, int_wavelet_array, 1);
	image_set_color_to_coordinates(image_scale, color, int_scale_array, 1);
	
	image_t * image_transformed = NULL;
	image_cat_hor(&image_transformed, image_scale, image_wavelet);
	image_save_ppm( image_transformed, "transformed.ppm");


	image_save_ppm(image, "image.ppm");
	
	image_free(image);

	return 0;
}

//Never tested the following function, should do some sort of circular convolution
//Keeping it there in case it comes handy in the future
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


