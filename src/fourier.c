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
#include <cyan/image/fourier.h>

//The function FFT_1D() is implemented as the procedure FFT(A) in the following link
//http://people.scs.carleton.ca/~maheshwa/courses/5703COMP/16Fall/FFT_Report.pdf
//FFT_1D computes (and returns in f) the 1D Fast Fourier Transform of the sequence f of 2^n elements (i.e. {f(k)}_{k=0,...,2^n}
//At this stage of development complexity is N * lgN in size and computations, it can be made in 2N size
//Data ordering is a little bit unclear
//At this point it is a big memory leak, will be fixed soon
complex_polar_t * FFT_1D( complex_polar_t * f, complex_polar_t * buffer, int n){
	
	int N = pow(2, n);
	if(buffer == NULL){
	//	fprintf(stdout, "FFT_1D : Allocating the buffer (of size N : %d).\n", N);
		buffer = (complex_polar_t *) malloc( N * sizeof(complex_polar_t));
	}
	if(n == 0)
		return f;
	
	complex_polar_t * W_N = unity(1, N);
	complex_polar_t * W = unity(0, 1);
	
	complex_polar_t * temp_polar = unity(0, 1);	
	complex_cart_t  temp_cart1;
	complex_cart_t  temp_cart2;

	complex_cart_t f_temp_cart;
	int i;
	//Rearranging f with the N/2 even elems, then the N/2 odd elems
	//(A_even)_{2k} corresponds to order_buffer[0, N/2-1]
	//(A_odd)_{2k+1} corresponds to order_buffer[N/2, N-1] 
	
	complex_polar_t order_buffer[N];
	for(i = 0; i< N/2; i++){
		order_buffer[i] = f[2*i]; //to N-2
		order_buffer[N/2+i] = f[2*i + 1];	
	}

	complex_polar_t * Y_even;
	complex_polar_t * Y_odd;
	Y_even = FFT_1D(order_buffer, buffer, n-1);
	Y_odd = FFT_1D(&(order_buffer[N/2]), &(buffer[N/2]), n-1);
	
	int j = 0;
	for(j = 0; j< N/2; j++)
	{
		//sorting values part
		mult_complex_polar(temp_polar, *W, Y_odd[ j]);	//temp_polar i.e W * Y_odd
		polar_to_cart(&temp_cart1, *temp_polar);	//temp_cart 1 i.e W*Y_odd
		polar_to_cart(&temp_cart2, Y_even[j]);	//temp_cart2 i.e. Y_even
		
		//Y(j) part

		add_complex_cart(&f_temp_cart, temp_cart1, temp_cart2);
		cart_to_polar(&(buffer[j]), f_temp_cart);
			
		//Y(j+N/2) part
		
		substract_complex_cart(&f_temp_cart, temp_cart2, temp_cart1);
		cart_to_polar(&(buffer[j + N/2]), f_temp_cart );
		
		//Updating phase factor
		//mult_complex_polar(temp_polar, *W, *W_N);
		W->phase += 2.0f * (4.0f * atan(1.0f)) / (double) N;	//Phase += 2*pi /N	
		//W->power = temp_polar->power;	
	}

	return buffer;
}


//Takes as an input fourier_polar, the polar fourier coefficients and computes f_polar the reverse FT of fourier_polar
//
complex_polar_t * FFT_1D_reverse( complex_polar_t * fourier_polar, int n){

//	rotate_buffer(fourier_polar, pow(2, n), sizeof(complex_polar_t));
	complex_polar_t * f_polar;
	f_polar = FFT_1D(fourier_polar, NULL, n);
	unitary_ft_polar( f_polar, pow(2, n));
	reverse_buffer(f_polar, pow(2, n), sizeof(complex_polar_t));

	return f_polar;
}

complex_cart_t * FFT_1D_cart_to_cart( complex_cart_t * f, int n){
	complex_polar_t * f_polar = NULL;
	int N = pow(2, n);
		
	f_polar = cart_array_to_polar(f, N);

	//rotate_buffer(f_polar, N, sizeof(complex_polar_t));	
	complex_polar_t * fft_polar = NULL;

	fft_polar = FFT_1D(f_polar, NULL, n);
	unitary_ft_polar( fft_polar, N);
//	rotate_buffer(fft_polar, N, sizeof(complex_polar_t));
	
	free(f_polar);

	complex_cart_t * fft_cart = NULL;
	fft_cart = polar_array_to_cart(fft_polar, N);

	free(fft_polar);

	return fft_cart;
}


complex_cart_t * FFT_1D_reverse_cart_to_cart( complex_cart_t * fft, int n){
	complex_polar_t * fft_polar = NULL;
	int N = pow(2, n);
		
	fft_polar = cart_array_to_polar(fft, N);
	
	complex_polar_t * f_polar = NULL;

//	rotate_buffer(fft_polar, N, sizeof(complex_polar_t));
	f_polar = FFT_1D_reverse(fft_polar, n);

	rotate_buffer(fft_polar, N, sizeof(complex_polar_t));
	//reverse_buffer(fft_polar, N, sizeof(complex_polar_t));
	free(fft_polar);

	complex_cart_t * f_cart = NULL;
	f_cart = polar_array_to_cart(f_polar, N);

	free(f_polar);

	return f_cart;

}

//Returns a double array from a 2^n lines by 2^m cols initial array
complex_cart_t ** FFT_2D(complex_cart_t ** array_cart, int n, int m){
	int N = pow(2, n);
	int M = pow(2, m);


	complex_cart_t * ft_array_rows[M];
	complex_cart_t temp_row[M];
	complex_cart_t ** ft_array_cols = NULL;
	ft_array_cols = (complex_cart_t **) malloc( N * sizeof( complex_cart_t * ) );
	int i, j, coords;
	for(j = 0; j < M; j++){	
		ft_array_rows[j] = FFT_1D_cart_to_cart(array_cart[j], m);
	}
	for(i=0; i < M; i++){
		for(j = 0; j < N; j ++){
			temp_row[j] = ft_array_rows[j][i];
		}	
		ft_array_cols[i] = FFT_1D_cart_to_cart( temp_row, n);
	}


	return ft_array_cols;
}
complex_cart_t ** FFT_2D_reverse(complex_cart_t *** dst, complex_cart_t ** ft_array_cart, int n, int m){
	int N = pow(2, n);
	int M = pow(2, m);

	complex_cart_t * ft_array_rows[N];
	complex_cart_t ft_temp_row[M];
	complex_cart_t ** array_cols = NULL;
	array_cols = (complex_cart_t **) malloc( M * sizeof( complex_cart_t * ) );
	int i, j, coords;
	for(i = 0; i < N; i++){
		ft_array_rows[i] = FFT_1D_reverse_cart_to_cart(ft_array_cart[i], m);
	}
	for(i = 0; i < N; i++){
		for(j=0; j < M; j++){
			ft_temp_row[j] = ft_array_rows[j][i];
		}
		array_cols[i] = FFT_1D_reverse_cart_to_cart(ft_temp_row, n);
	}	

//	rotate_180( (void ***) dst, (void **)  array_cols, N, sizeof(complex_cart_t * ) );

	return *dst = array_cols;
}

//unitary_ft_*() divide each coefficients inside *_array by 1/sqrt(N) in order to have a unitary fourier transform.
void * unitary_ft_polar(complex_polar_t * polar_array, int N){
	int i;
	for(i=0; i<N; i++)
		polar_array[i].power *= 1.0f/sqrt( (double) N);
}

void * unitary_ft_cart(complex_cart_t * cart_array, int N){
	int i;
	for(i=0; i<N; i++){
		cart_array[i].real *= 1.0f/sqrt( (double) N);
		cart_array[i].im *= 1.0f/sqrt( (double) N);
	}
}




