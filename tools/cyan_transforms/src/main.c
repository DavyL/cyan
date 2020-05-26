#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#include <cyan/color/color.h>
#include <cyan/image/image.h>

#include <cyan_fileio/load_png.h>
#include <cyan_fileio/save_ppm.h>

#include <cyan/image/complex.h>
#include <cyan/image/fourier.h>
#include <cyan/image/transforms.h>
#include <cyan/image/sampling.h>

int main( int argc, char** argv, char* envv ) {

	int i, j ;
	int result ;
	image_t * image = NULL ;	
	
	if (argc != 2 ) {
		fprintf(stderr,"Usage : %s file.png\n", argv[0] ) ;
		return -1 ;
	}
	FILE * fp;

	fp = fopen(argv[1], "r");
	if(fp == NULL){	
		fprintf(stderr, "Couldn't open file.\n Error : %d, (%s)\n", errno, strerror(errno));
		return -1;
	}
	image = png2image(fp);
//	fclose(fp);

	if ( image == NULL ) {
		fprintf(stderr,"image allocation failed\n") ;
		return -1 ;
	}	


	/////////Computing fourier transform

	int n, m;

	n = (int) log2( (double) image->rows );
	m = (int) log2( (double) image->cols );
	
	if( (int) pow(2, n + m) != image->rows * image->cols ){
		fprintf(stderr, "ERROR: main: image : %s \t is not a product of powers of two.\n");
	       return -1;	
	}


	//Converting the image in a complex array
	//Y component is sent to the real part
	complex_cart_t ** image_array = NULL;
	image_array = image_to_cart(image, Y_to_cart); 


	complex_cart_t ** window_array = NULL;

	array_2d_through_window_arb( (void ***) &window_array, (void **) image_array, image->rows, image->cols, sizeof(complex_cart_t),
			window_triangular, mult_scalar_complex_cart);	


	image_t * grey_image_window = NULL;
	grey_image_window = cart_to_Y(window_array, image->rows, image->cols);
	image_save_ppm(grey_image_window, "window_grey_image.ppm");

	image_t * grey_image = NULL;
	grey_image = cart_to_Y(image_array, image->rows, image->cols);
	image_save_ppm(grey_image, "grey_image.ppm");

	//Computing Fourier transform
	complex_cart_t ** ft_array = NULL;
	ft_array = FFT_2D(image_array, n, m);

	
	///////Saving various images
	
	
	////Saving image of fourier transform
	image_t * ft_image = NULL;
	ft_image = cart_to_Y(ft_array, image->rows, image->cols);	
	image_save_ppm ( ft_image, "fourier.ppm");

	//Saves power coefficients of fourier coeficients in an image
	image_t * power_image = NULL;
	power_image = cart_power_to_Y(ft_array, image->rows, image->cols );
	image_save_ppm( power_image, "power.ppm" );
	image_free( power_image);

	//Saves the log of the power of fourier coeficients in an image
	image_t * log_power_image = NULL;
	log_power_image = cart_log_power_to_Y(ft_array, image->rows, image->cols);
	image_save_ppm( log_power_image, "log_power.ppm" );
	image_free( log_power_image );

	//Saves the phase of fourier coeficients in an image
	image_t * phase_image = NULL;
	phase_image = cart_phase_to_Y( ft_array, image->rows, image->cols );
	image_save_ppm( phase_image, "phase.ppm" );
	image_free( phase_image );

	////////Reverting the fourier coeficients

	//Computes the reverse FT on the coeficients and saves it in an image
	image_t * reverse_image = NULL;
	//reverse_image = FT_image_Y( ft_image, FFT_2D_reverse, cart_to_Y);
	
	complex_cart_t ** reverse_array = NULL	;
//	rotate_buffer( ft_array, image->rows, sizeof(complex_cart_t *));

	FFT_2D_reverse(&reverse_array, ft_array, n, m);

	reverse_image = cart_to_Y ( reverse_array, image->rows, image->cols );
	image_save_ppm(reverse_image, "reverse_fourier.ppm");
	image_free(reverse_image);

	//image_free( image );

	




	
	image_free(image);
	return 0 ;
}


