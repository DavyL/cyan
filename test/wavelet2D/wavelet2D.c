#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cyan/image/image.h>
#include <cyan/image/fourier.h>
#include <cyan/image/sampling.h>
#include <cyan/image/transforms.h>
#include <cyan/image/wavelet.h>

#include <cyan/color/color.h>

#include <cyan_fileio/load_png.h>
#include <cyan_fileio/save_ppm.h>

int transform_2_dim( double ** dst1, double ** dst2, double ** dst3, double ** dst4, double * src, int len, int height, int (*transform)(double **, double **, double *, int , int));
int transform_2_dim_backward( double ** dst1, double * src1, double * src2, double * src3, double * src4, int len, int height, int (*transform)(double **, double *, double *, int , int));

int save_wavelet_image(image_t ** dst, 	double * sca_sca, double * wav_sca, double * sca_wav, double * wav_wav, int rows, int cols);

int main( int argc, char** argv, char* envv ) {
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

	int len = image->cols;
	int height = image->rows;
	int offset = 1;
	int i = 0;
	int j = 0;
	int scale = 255;
	double array[len];

	double * sca_sca = NULL;
	double * wav_sca = NULL;
	double * sca_wav = NULL;
	double * wav_wav = NULL;

	D4_wavelet_2D_CR( &sca_sca, &sca_wav, &wav_sca, &wav_wav, image->X, image->rows, image->cols);

	image_t * image_transformed = NULL;
	save_wavelet_image(&image_transformed, sca_sca, wav_sca, sca_wav, wav_wav, image->rows, image->cols);

	
/*
	//Non-linear approx
	double T = 0.2f;
	double_threshold(scale_array, T, len/2);
	double_threshold(wavelet_array, T, len/2);
*/

	double array_norm = 0.0;
	l2_distance(&array_norm, image->X, NULL, image->cols * image->rows );
	fprintf(stdout, " array norm : %f\n", array_norm );
	double norm[4];
	l2_distance(norm, sca_sca, NULL, image->cols * image -> rows /4); 
	l2_distance(norm + 1, sca_wav, NULL, image->cols * image -> rows /4); 
	l2_distance(norm + 2, wav_sca, NULL, image->cols * image -> rows /4); 
	l2_distance(norm + 3, wav_wav, NULL, image->cols * image -> rows /4); 
	
	fprintf(stdout, "norms : %f\t%f\t%f\t%f\n", norm[0], norm[1], norm[2], norm[3]);	
	
	image_save_ppm(image, "original.ppm");
	image_save_ppm(image_transformed, "transformed.ppm");	



	return 0;
}
int save_wavelet_image(image_t ** dst, 	double * sca_sca, double * wav_sca, double * sca_wav, double * wav_wav, int rows, int cols){

	if(dst == NULL || sca_sca == NULL || wav_sca == NULL || sca_wav == NULL || wav_wav == NULL){
		fprintf(stdout, "ERR : save_wavelet_image() : One of the arguments is a NULL pointer.\n");
		return -1;
	}
	
	image_t * image_sca_sca = image_new(rows/2, cols/2);
	image_sca_sca->X = sca_sca;
	image_t * image_wav_sca = image_new(rows/2, cols/2);
	image_wav_sca->X = wav_sca;
	image_t * image_scale = NULL;
	image_cat_hor(&image_scale, image_wav_sca, image_sca_sca);
	
	image_t * image_wav_wav = image_new(rows/2, cols/2);
	image_wav_wav->X = wav_wav;
	image_t * image_sca_wav = image_new(rows/2, cols/2);
	image_sca_wav->X = sca_wav;
	image_t * image_wavelet = NULL;	
	image_cat_hor(&image_wavelet, image_wav_wav, image_sca_wav);

	image_cat_ver(dst, image_scale, image_wavelet);
/*	Freeing the images also frees the wavelet input
	image_free(image_sca_sca);	
	image_free(image_sca_wav);	
	image_free(image_wav_sca);	
	image_free(image_wav_wav);
*/
	image_free(image_scale);
	image_free(image_wavelet);

	return 0;
	

}
//Not sure this is the right way to do it
//Correspondance of dst* :
//1 : psi_y_psi_x 
//2 : phi_y_psi_x 
//3 : psi_y_phi_x 
//4 : phi_y_phi_x 
//psi: lowpass
//phi: highpass (details)
int transform_2_dim( double ** dst1, double ** dst2, double **dst3, double ** dst4, double * src, int src_len, int src_height, int (*transform)(double **, double **, double *, int , int)){

	int i = 0;
	int len = src_len / 2;
	double * xwavelet = (double *) malloc(len * src_height * sizeof(double));
	double * xscale = (double *) malloc( len * src_height * sizeof(double));

	for( i = 0; i < src_height* src_len; i++){
		fprintf(stdout, "%f\t", *(src + i));
	}
	double *p1, *p2, * p3;
	for( i = 0; i < src_height; i++){
		p1 = xwavelet + i * len ;
		p2 = xscale + i * len ;
		p3 = src + i * src_len ;

		transform(&p1, &p2, p3, src_len, 1);
		//fprintf(stdout, "p1 : %f\t p2 : %f\t p3%f\n", *p1, *p2, *p3);
//		fprintf(stdout, "i : %d\nwav: %f\t%f\t%f\nsca: %f\t%f\t%f\n \n", i, *p1, p1[10], p1[20], p2[0], p2[10], p2[20]);	
	}
	
	int height = src_height / 2;
	if(*dst3 == NULL){
		*dst3 = (double *) malloc(len * height * sizeof(double));
	}
	if(*dst4 == NULL){
		*dst4 = (double *) malloc(len * height * sizeof(double));
	}

	for( i = 0; i < len; i++){
		p1 = *dst3 + i ;
		p2 = *dst4 + i ;
		p3 = xwavelet + i ;

		transform(&p1, &p2, p3, src_height , len);
		//fprintf(stdout, "p1 : %f\t p2 : %f\t p3%f\n", *p1, *p2, *p3);
	}
	//*dst3 = psi_y_phi_x;
	//*dst4 = phi_y_phi_x;
	if(*dst1 == NULL){
		*dst1 = (double *) malloc(len * height * sizeof(double));
	}
	if(*dst2 == NULL){
		*dst2 = (double *) malloc(len * height * sizeof(double));
	}

	for( i = 0; i < len; i++){
		p1 = *dst1 + i;
		p2 = *dst2 + i;
		p3 = xscale + i;

		transform(&p1, &p2, p3, src_height, len);
		//fprintf(stdout, "p1 : %f\t p2 : %f\t p3%f\n", *p1, *p2, *p3);
	}
	//*dst1 = psi_y_psi_x;
	//*dst2 = phi_y_psi_x;

	free(xwavelet);
	free(xscale);

	return 0;
}

//See the forward version for correspondance of src*
//len and height correspond to the size of src*
//*dst1 is of size 4*|src*|
int transform_2_dim_backward( double ** dst1, double * src1, double * src2, double * src3, double * src4, int len, int height, int (*transform)(double **, double *, double *, int , int)){
		
	int i = 0;
	double * p1 = NULL;	//Temporary pointers to pass as arguments
	double * p2 = NULL;
	double * p3 = NULL;

	//Compute x analysis of the image
	double * scale_x = malloc(len * 2 * height * sizeof(double) );
	for(i = 0; i < len; i++){
		p1 = scale_x + i;
		p2 = src3 + i;
		p3 = src4 + i;
		transform(&p1, p2, p3, height, len);
		//fprintf(stdout, "p1 : %f\t p2 : %f\t p3%f\n", *p1, *p2, *p3);
	}

	double * wavelet_x = malloc(len * 2 * height * sizeof(double) );

	for(i = 0; i < len; i++){
		p1 = wavelet_x + i;
		p2 = src1 + i;
		p3 = src2 + i;
		transform(&p1, p2, p3, height, len);

	}

	//Merging the two x analysis to produce the original image
	
	if(*dst1 == NULL){
		*dst1 = malloc( 2 * len * 2 * height * sizeof(double));
	}
	for(i = 0; i < 2 * height; i++){
		p1 = *dst1 + 2 * i * len;
		p2 = wavelet_x + i * len;
		p3 = scale_x + i * len;
		
		transform(&p1, p2, p3, len, 1);
	}

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


