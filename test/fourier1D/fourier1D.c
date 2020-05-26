#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#include <cyan/color/color.h>
#include <cyan/image/image.h>
#include <cyan/image/load_png.h>
#include <cyan/image/transforms.h>
#include <cyan/image/fourier.h>

//This unit test  computes the fourier transform of cart_array
//Then the reverse ft of the array of coefs. obtained
//The values are taken from the example that can be found on the wikipedia page of the DFT

int main( int argc, char** argv, char* envv ) {

	int i, j ;
	complex_cart_t cart_array[4];

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
	reverse_buffer(reverse_cart, N, sizeof(complex_cart_t));

	for(i = 0; i < N; i++){
		fprintf(stdout, "array reverse : real : %f, im : %f \n", reverse_cart[i].real, reverse_cart[i].im);
	}

	return 0 ;
}


