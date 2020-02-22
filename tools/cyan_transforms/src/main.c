#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#include <cyan/color/color.h>
#include <cyan/image/image.h>
#include <cyan/image/load_png.h>

#include "transforms.h"

int main( int argc, char** argv, char* envv ) {

	int i, j ;
	int result ;
	image_t * image ;	
	image_t * grey_image;
	FILE * fp;

	fp = fopen("lena.png", "r");
	if(fp == NULL){	
		fprintf(stderr, "Couldn't open file.\n Error : %d, (%s)\n", errno, strerror(errno));
		return -1;
	}


	if (argc != 2 ) {
		fprintf(stderr,"Usage : %s output_file.ppm\n", argv[0] ) ;
		return -1 ;
	}
	image = png2image(fp);
	if ( image == NULL ) {
		fprintf(stderr,"image allocation failed\n") ;
		return -1 ;
	}	
//	grey_image = color2grey(image);

//	result = image_save_ppm( grey_image, argv[1] ) ; 

	complex_cart_t cart_array[4];
	//complex_polar_t array[4];
	//complex_polar_t buffer[4];

	cart_array[0].real = 1.0f;
	cart_array[0].im = 0.0f;
	
	cart_array[1].real = 2.0f;
	cart_array[1].im = -1.0f;
	
	cart_array[2].real = 0.0f;
	cart_array[2].im = -1.0f;
	
	cart_array[3].real = -1.0f;
	cart_array[3].im = 2.0f;

	int n = 2;
	int N = pow(2, n);
	complex_cart_t * fft_cart = NULL;
	fft_cart = FFT_1D_cart_to_cart(cart_array, n);
	for(i=0; i < N; i++){
		fprintf(stdout, "array fft : real : %f, im : %f \n", fft_cart[i].real, fft_cart[i].im);
	}

	complex_cart_t * reverse_cart = NULL;
	reverse_cart = FFT_1D_reverse_cart_to_cart(fft_cart, n);
	for(i = 0; i < N; i++){
		fprintf(stdout, "array reverse : real : %f, im : %f \n", reverse_cart[i].real, reverse_cart[i].im);
	}

	if (result != 0) {
		fprintf(stderr, "cannot save file\n") ;
		return -1 ;
	}

	image_free( image) ;

	fclose(fp);

	return 0 ;
}


