#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cyan/image/image.h>
#include <cyan/image/fourier.h>
#include <cyan/image/sampling.h>

#include <cyan/color/color.h>

#include <cyan_fileio/save_ppm.h>

int main( int argc, char** argv, char* envv ) {

	int len = 256;
	int i = 0;
	int scale = 255;
	double array[len];
	complex_polar_t polar_array[len];

	for(i=0; i < len; i++){
		//array[i] = cos( 2 * M_PI * (double) i / (double) len ) + 
		array[i] =  3.0f * (double) i / (double) len;
		polar_array[i].phase = 0.0f;
		polar_array[i].power = array[i];
	}

	complex_polar_t * ft = NULL;	
	ft = FFT_1D( polar_array, NULL, 8);
	unitary_ft_polar( ft, len);
	
	
	double phase_array[len];
	double power_array[len];
	double rev_phase_array[len];
	double rev_power_array[len];

	complex_polar_t * rev_ft = NULL;
	rev_ft = FFT_1D_reverse(ft, 8);

	rotate_buffer( ft, len, sizeof(complex_polar_t));



	for(i=0; i < len; i++){
		phase_array[i] = ft[i].phase;
		power_array[i] = ft[i].power;
		rev_phase_array[i] = rev_ft[i].phase;
		rev_power_array[i] = rev_ft[i].power;
		fprintf(stdout, " rev_ft[%d].power = %f\n", i, rev_ft[i].power);
	}


	normalize_and_scale_double_array(array, len, (double) scale);
	normalize_and_scale_double_array(phase_array, len, (double) scale);
	normalize_and_scale_double_array(power_array, len, (double) scale);
	normalize_and_scale_double_array(rev_phase_array, len, (double) scale);
	normalize_and_scale_double_array(rev_power_array, len, (double) scale);



	int * int_array = NULL;
	int * int_phase_array = NULL;
	int * int_power_array = NULL;
	double_array_to_int_array( &int_phase_array, phase_array, len);
	double_array_to_int_array( &int_power_array, power_array, len);
	
	int * int_rev_phase_array = NULL;
	int * int_rev_power_array = NULL;
	double_array_to_int_array( &int_rev_phase_array, rev_phase_array, len);
	double_array_to_int_array( &int_rev_power_array, rev_power_array, len);
	
	double_array_to_int_array( &int_array, array, len);

	image_t * image = image_new_empty(len, scale);
	
	color_t color;
	color_assign( &color, 1.0f, 1.0f, 1.0f, 0);

	image_set_color_to_coordinates(image, color, int_array, 1);
	
	color_assign( &color, 1.0f, 1.0f, 1.0f, 0);
	//image_set_color_to_coordinates(image, color, int_phase_array, 1);
	
	color_assign( &color, 1.0f, 0.5f, 0.5f, 0);
	image_set_color_to_coordinates(image, color, int_power_array, 1);
	image_set_color_to_coordinates(image, color, int_rev_power_array, 1);
	
	image_save_ppm(image, "image.ppm");
	
//	image_free(image);

	return 0;

}
