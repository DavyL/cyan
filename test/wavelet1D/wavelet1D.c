#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cyan/image/image.h>
#include <cyan/image/fourier.h>
#include <cyan/image/sampling.h>
#include <cyan/image/transforms.h>

#include <cyan/color/color.h>

#include <cyan_fileio/save_ppm.h>

int convol_loc( double * a, double * array, double * h, int i, int filter_len, int len );

int main( int argc, char** argv, char* envv ) {

	int len = 256;
	int i = 0;
	int scale = 255;
	double array[len];
	complex_polar_t polar_array[len];
	
	//Loading Daubechies4 low pass filter
	const int filter_length = 4;
	double h[] = { 0.482962913f, 0.836516303f, 0.224143868f, -0.129409522f}; 
	//Computing the associated high pass filter
	//TODO : implement general computation of g from h
	double g[] = { -0.12940952f, -0.22414387f, 0.8365163f, -0.48296291f};

	mult_double_array(h, filter_length, 1.0f);
	mult_double_array(g, filter_length, 1.0f);
	for(i=0; i < len; i++){
		//array[i] = cos( 2 * M_PI * (double) i / (double) len ) + 
		array[i] =  3.0f * (double) i / (double) len;
		polar_array[i].phase = 0.0f;
		polar_array[i].power = array[i];
	}
	
	double a[len];
	double d[len];
	double temp = 0.0f;
	for(i = 0; i < len; i++){
		convol_loc( &temp, array, h, i, filter_length, len ); 
		a[i] = temp;
		fprintf(stdout, "a %d : %f\n", i, temp);
		convol_loc( &temp, array, g, i, filter_length, len ); 
		d[i] = temp;
	}
	double padded_a[len/2];
	double padded_d[len/2];
	for(i=0; i < len/2; i++){
		padded_a[i] = a[2*i];
		//padded_a[2*i + 1] = 0.0f;
	
		padded_d[i] = d[2*i];
		//padded_d[2*i + 1] = 0.0f;
	}
	//reverse_buffer(h, filter_length, sizeof(double));
	//reverse_buffer(g, filter_length, sizeof(double));
	double reverse_h[filter_length]; 
	reverse_h[0] = h[2];
	reverse_h[1] = g[2]; 
	reverse_h[2] = h[0];
	reverse_h[3] = g[0];
	double reverse_g[filter_length];
	reverse_g[0] = h[3];
	reverse_g[1] = g[3]; 
	reverse_g[2] = h[1];
	reverse_g[3] = g[1];
	
	double reversed_array[len];
	double detail_temp = 0;
	double harsh_temp = 0;
	int j = 0;
	reversed_array[j++] = a[len - 1] * reverse_h[0] + d[len - 1] * reverse_h[1] + a[0]*reverse_h[2] + d[0]*reverse_h[3];
	reversed_array[j++] = a[len - 1] * reverse_g[0] + d[len - 1] * reverse_g[1] + a[0]*reverse_g[2] + d[0]*reverse_g[3];
	for(i = 0; i < len/2 ; i++){
		reversed_array[j++] = a[2*i] * reverse_h[0] + d[2*i] * reverse_h[1] + a[2*i+2]*reverse_h[2] + d[2*i + 2]*reverse_h[3];
		reversed_array[j++] = a[2*i] * reverse_g[0] + d[2*i] * reverse_g[1] + a[2*i+2]*reverse_g[2] + d[2*i + 2]*reverse_g[3];

	}

	double L2_diff = 0;
	for( i = 0; i < len; i++){
		fprintf(stdout, " array : %f, \t reversed array : %f\n", array[i], reversed_array[i]);
		L2_diff += fabs( array[i]*array[i] - reversed_array[i]*reversed_array[i] );
	}
	fprintf(stdout, "L2 difference of array and reversed array is : %f.\n", L2_diff);

	int * int_array = NULL;
	int * int_reversed_array = NULL;
	int * int_a	= NULL;
	int * int_d	= NULL;

	normalize_and_scale_double_array(array, len, (double) scale);
	normalize_and_scale_double_array(reversed_array, len, (double) scale);
	normalize_and_scale_double_array(a, len/2, (double) scale);
	normalize_and_scale_double_array(d, len/2, (double) scale);
	double_array_to_int_array( &int_array, array, len);
	double_array_to_int_array( &int_reversed_array, reversed_array, len);
	double_array_to_int_array( &int_a, a, len/2);
	double_array_to_int_array( &int_d, d, len/2);

	image_t * image = image_new_empty(len, scale);
	
	color_t color;
	color_assign( &color, 1.0f, 1.0f, 1.0f, 0);

	image_set_color_to_coordinates(image, color, int_array, 1);
	
	color_assign( &color, 0.7f, 0.3f, 0.8f, 0);

	image_set_color_to_coordinates(image, color, int_reversed_array, 1);
	
	color_assign( &color, 0.5f, 1.0f, 1.0f, 0);
	//image_set_color_to_coordinates(image, color, int_a, 1);
	
	color_assign( &color, 1.0f, 1.0f, 0.5f, 0);
	//image_set_color_to_coordinates(image, color, int_d, 1);


	image_save_ppm(image, "image.ppm");
	
//	image_free(image);

	return 0;
}
		
//Stores in a the convolution of array (of length n) 
//at step i with a filter h of length filter_len
int convol_loc( double * a, double * array, double * h, int i, int filter_len, int len ){
	if( a == NULL || array == NULL || h == NULL){
		fprintf(stderr, "ERR : wavelet1D : convol() : a NULL pointer was passed as an argument.\n");
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


