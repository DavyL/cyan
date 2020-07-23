#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <cyan/image/transforms.h>
#include <cyan/image/wavelet.h>

//Computes the Daubechies 4 forward wavelet transform with as in input, src, a double array of length len
//Output is stored in dst_*, which are arrays of length len / 2
int D4_wavelet_forward( double ** dst_scale, double ** dst_wavelet, double * src, int len, int offset){

	if(dst_scale == NULL || dst_wavelet == NULL || src == NULL){
		fprintf(stderr, "Error : D4_wavelet_forward() A NULL argument was passed.\n");
		return -1;
	}
	if(len % 2){
		fprintf(stderr, "Error : D4_wavelet_forward() len is an odd number.\n");
		return -1;
	}
	if(*dst_scale == NULL){
		fprintf(stdout, "D4_wavelet_forward() : Allocating memory to dst_scale.\n");
		*dst_scale = malloc((len/2) * sizeof(double));
	}
	if(*dst_wavelet == NULL){
		fprintf(stdout, "D4_wavelet_forward() : Allocating memory to dst_wavelet.\n");
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

	int i = 0;	
	
	double temp = 0.0f;
	//Doing finite support convolutions w.r.t the associated filter
	for(i = 0; i < len / 2 ; i++){	
		convol_loc_1D( (*dst_scale) + i, src, h, 2*i, filter_length, len ); 
		convol_loc_1D( (*dst_wavelet) + i, src, g, 2*i, filter_length, len ); 
	}

	return 0;
}
//Computes the daubechies 4 wavelet transform of an image with len_row rows, len_col columns
//Hence src is an array of len_row * len_col elements
//The output is stored in dst_* of size (len_row/2) * len_col
int D4_wavelet_columns( double ** dst_scale, double ** dst_wavelet, double * src, int len_row, int len_col){

	if(dst_scale == NULL || dst_wavelet == NULL || src == NULL){
		fprintf(stderr, "Error : D4_wavelet_forward_columns() A NULL argument was passed.\n");
		return -1;
	}
	if(len_row % 2){
		fprintf(stderr, "Error : D4_wavelet_forward_columns() len is an odd number.\n");
		return -1;
	}
	if(*dst_scale == NULL){
		fprintf(stdout, "D4_wavelet_forward_columns() : Allocating memory to dst_scale.\n");
		*dst_scale = malloc((len_row/2)* len_col * sizeof(double));
	}
	if(*dst_wavelet == NULL){
		fprintf(stdout, "D4_wavelet_forward_columns() : Allocating memory to dst_wavelet.\n");
		*dst_wavelet = malloc((len_row/2)* len_col * sizeof(double));
	}
	if(*dst_wavelet == NULL || *dst_scale == NULL){
		fprintf(stderr, "Error: D4_wavelet_forward_columns() Error allocating dst pointers.\n");
		return -1;
	}

	//Loading Daubechies4 low pass filter (scaling function)
	const int filter_length = 4;
	double h[] = { 0.482962913f, 0.836516303f, 0.224143868f, -0.129409522f}; 
	
	//Loading the associated high pass filter (wavelet function)
	//TODO : implement general computation of g from h
	double g[] = { -0.12940952f, -0.22414387f, 0.8365163f, -0.48296291f};

	int i = 0;	
	int j = 0;

	double temp = 0.0f;
	//Doing finite support convolutions w.r.t the associated filter
	for(i = 0; i < len_row / 2; i++){
		for(j = 0; j < len_col; j++){	
			convol_loc_1D_columns( (*dst_scale) + i*len_col + j, src, h, 2*i, j, filter_length, len_row, len_col ); 
			convol_loc_1D_columns( (*dst_wavelet) + i*len_col + j, src, g, 2*i, j, filter_length, len_row, len_col );  
		}
	}

	return 0;
}
//Computes the daubechies 4 wavelet transform on the rows of an image with len_row rows, len_col columns
//Hence src is an array of len_row * len_col elements
//The output is stored in dst_* of size len_row * (len_col/2)
int D4_wavelet_rows( double ** dst_scale, double ** dst_wavelet, double * src, int len_row, int len_col){

	if(dst_scale == NULL || dst_wavelet == NULL || src == NULL){
		fprintf(stderr, "Error : D4_wavelet_forward_rows() A NULL argument was passed.\n");
		return -1;
	}
	if(len_row % 2){
		fprintf(stderr, "Error : D4_wavelet_forward_rows() len is an odd number.\n");
		return -1;
	}
	if(*dst_scale == NULL){
		fprintf(stdout, "D4_wavelet_forward_rows() : Allocating memory to dst_scale.\n");
		*dst_scale = malloc((len_row/2)* len_col * sizeof(double));
	}
	if(*dst_wavelet == NULL){
		fprintf(stdout, "D4_wavelet_forward_rows() : Allocating memory to dst_wavelet.\n");
		*dst_wavelet = malloc((len_row/2)* len_col * sizeof(double));
	}
	if(*dst_wavelet == NULL || *dst_scale == NULL){
		fprintf(stderr, "Error: D4_wavelet_forward_rows() Error allocating dst pointers.\n");
		return -1;
	}

	//Loading Daubechies4 low pass filter (scaling function)
	const int filter_length = 4;
	double h[] = { 0.482962913f, 0.836516303f, 0.224143868f, -0.129409522f}; 
	
	//Loading the associated high pass filter (wavelet function)
	//TODO : implement general computation of g from h
	double g[] = { -0.12940952f, -0.22414387f, 0.8365163f, -0.48296291f};


	int i = 0;	
	int j = 0;	
	double temp = 0.0f;
	//Doing finite support convolutions w.r.t the associated filter
	for(i = 0; i < len_row ; i++){
		for(j = 0; j < len_col/2; j++){	
			convol_loc_1D( (*dst_scale) + i * (len_col/2) + j, src + i*len_col, h, 2 * j, filter_length, len_col ); 
			convol_loc_1D( (*dst_wavelet) + i *(len_col/2) + j, src + i*len_col, g,  2 * j, filter_length, len_col );  
		}
	}

	return 0;
}


//Performs 2D wavelet (forward) transform of array of size len_row*len_col
//Returns wavelet-scale coefficient arrays of size (len_row/2) * (len_col/2)
int D4_wavelet_2D( double ** sca_sca, double ** sca_wav, double ** wav_sca, double ** wav_wav, double * array, int len_row, int len_col){


	
	double * hor_scale = NULL;
	double * hor_wavelet = NULL;

	D4_wavelet_rows( &hor_scale, &hor_wavelet, array, len_row, len_col);
	D4_wavelet_columns( sca_sca, sca_wav, hor_scale, len_row, len_col/2);
	D4_wavelet_columns( wav_sca, wav_wav, hor_wavelet, len_row, len_col/2);


	return 0;
}

//Same as D4_wavelet_2D but performs Columns first, then Rows
int D4_wavelet_2D_CR( double ** sca_sca, double ** sca_wav, double ** wav_sca, double ** wav_wav, double * array, int len_row, int len_col){
	
	double * ver_scale = NULL;
	double * ver_wavelet = NULL;

	D4_wavelet_columns( &ver_scale, &ver_wavelet, array, len_row, len_col);
	D4_wavelet_rows( sca_sca, sca_wav, ver_scale, len_row/2, len_col);
	D4_wavelet_rows( wav_sca, wav_wav, ver_wavelet, len_row/2, len_col);
	return 0;
}
//Takes as an input scale and wavelet arrays given by forward Daubechies 4 wavelet transform of length len
//stores in dst the resulting backward transform of length 2 * len	
int D4_wavelet_backward( double ** dst, double * src_scale, double * src_wavelet, int len, int offset){
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
			fprintf(stderr, "Error: D4_wavelet_backward() : Couldn't allocate memory.\n");
			return -1;
		}
	}

	//reverse_* are reversed (fourier) coefs from * (h or g)
	double reverse_h[] = {-0.129409522, 0.224143868, 0.83651630, 0.482962913};
	double reverse_g[] = {-0.48296291, 0.8365163, -0.22414387, -0.12940952};

	//To recover signal we must compute scale_{j+1}  = (scale ^2 ) * h + (wavelet ^2) * g
	//(With ^2 the upsampling (by inserting 0s) and * the convolution 
	double conj_h[] = { reverse_h[0], reverse_h[2], reverse_h[1], reverse_h[3] };
	double conj_g[] = { reverse_g[0], reverse_g[2], reverse_g[1], reverse_g[3] };

	int i = 0;
	double tmp = 0.0;
	for(i = 0; i < len; i++){	
		convol_loc_1D( (*dst) + 2*i, src_scale, conj_h, i, 2, len ) + convol_loc_1D((*dst) + 2*i, src_wavelet, conj_g, i, 2, len); 
		convol_loc_1D( (*dst) + 2*i + 1, src_scale, &(conj_h[2]), i, 2, len ) + convol_loc_1D((*dst) + 2*i + 1, src_wavelet, &(conj_g[2]), i, 2, len); 
	}

	return 0;
}
//Takes as an input scale and wavelet (2D) arrays  given by forward Daubechies 4 wavelet transform of length len_row*len_col
//stores in dst the resulting backward transform on rows of length len_row * (2 * len_col)	
int D4_wavelet_backward_rows( double ** dst, double * src_scale, double * src_wavelet, int len_row, int len_col){
	if(dst == NULL || src_scale == NULL || src_wavelet == NULL){
		fprintf(stderr, "Invalid argument : D4_wavelet_backward() One of the arguments is a NULL pointer.\n");
		return -1;
	}
	if(len_row < 2){
		fprintf(stderr, "Invalid argument : D4_wavelet_backward() len is than 2.\n");
		return -1;
	}
	if(*dst == NULL){
		*dst = malloc(len_row * ( 2 * len_col ) * sizeof(double));
		if(*dst == NULL){
			fprintf(stderr, "Error: D4_wavelet_backward() : Couldn't allocate memory.\n");
			return -1;
		}
	}
	//Loading Daubechies4 low pass filter (scaling function)
	const int filter_length = 4;
	double h[] = { 0.482962913f, 0.836516303f, 0.224143868f, -0.129409522f}; 
	
	//Loading the associated high pass filter (wavelet function)
	//TODO : implement general computation of g from h
	double g[] = { -0.12940952f, -0.22414387f, 0.8365163f, -0.48296291f};



	double * dual_h = NULL;
	double * dual_g = NULL;

	compute_dual_filter(&dual_h, h, 4);
	compute_dual_filter(&dual_g, g, 4);

	//To recover signal we must compute scale_{j+1}  = (scale ^2 ) * h + (wavelet ^2) * g
	//(With ^2 the upsampling (by inserting 0s) and * the convolution 
	int i = 0;
	int j = 0;
	double tmp = 0.0;
	for(i = 0; i < len_row; i++){	
		for(j = 0; j < len_col; j++){
			convol_loc_1D( (*dst) + i* 2 * len_col + 2*j, src_scale + i * len_col, dual_h, j, 2, len_col); 
			convol_loc_1D((*dst) + i* 2 * len_col + 2*j, src_wavelet + i*len_col, dual_g, j, 2, len_col); 

			convol_loc_1D( (*dst) + i* 2 * len_col + 2*j +1, src_scale + i * len_col, &(dual_h[2]), j, 2, len_col );
			convol_loc_1D((*dst) + i * 2 * len_col + 2*j +1, src_wavelet + i*len_col, &(dual_g[2]), j, 2, len_col); 
		}
	}

	return 0;
}

//Takes as an input scale and wavelet (2D) arrays  given by forward Daubechies 4 wavelet transform of length len_row*len_col
//stores in dst the resulting backward transform on columns of length (2*len_row)*len_col	
int D4_wavelet_backward_columns( double ** dst, double * src_scale, double * src_wavelet, int len_row, int len_col){
	if(dst == NULL || src_scale == NULL || src_wavelet == NULL){
		fprintf(stderr, "Invalid argument : D4_wavelet_backward() One of the arguments is a NULL pointer.\n");
		return -1;
	}
	if(len_row < 2){
		fprintf(stderr, "Invalid argument : D4_wavelet_backward() len is than 2.\n");
		return -1;
	}
	if(*dst == NULL){
		*dst = malloc((2*len_row) * len_col * sizeof(double));
		if(*dst == NULL){
			fprintf(stderr, "Error: D4_wavelet_backward() : Couldn't allocate memory.\n");
			return -1;
		}
	}
	//Loading Daubechies4 low pass filter (scaling function)
	const int filter_length = 4;
	double h[] = { 0.482962913f, 0.836516303f, 0.224143868f, -0.129409522f}; 
	
	//Loading the associated high pass filter (wavelet function)
	//TODO : implement general computation of g from h
	double g[] = { -0.12940952f, -0.22414387f, 0.8365163f, -0.48296291f};



	double * dual_h = NULL;
	double * dual_g = NULL;

	compute_dual_filter(&dual_h, h, 4);
	compute_dual_filter(&dual_g, g, 4);

	//To recover signal we must compute scale_{j+1}  = (scale ^2 ) * h + (wavelet ^2) * g
	//(With ^2 the upsampling (by inserting 0s) and * the convolution 
	int i = 0;
	int j = 0;
	 
	double tmp = 0.0;
	for(i = 0; i < len_row; i++){	
		for(j = 0; j < len_col; j++){
			convol_loc_1D_columns( (*dst) + 2*i*len_col + j, src_scale, dual_h, i, j, 2, len_row, len_col); 
			convol_loc_1D_columns((*dst) + 2*i*len_col + j, src_wavelet, dual_g, i, j, 2, len_row, len_col); 

			convol_loc_1D_columns( (*dst) + (2*i + 1)*len_col + j, src_scale, &(dual_h[2]), i, j, 2, len_row, len_col ); 
			convol_loc_1D_columns((*dst) + (2*i + 1)*len_col + j, src_wavelet, &(dual_g[2]), i, j, 2, len_row, len_col); 
		}
	}

	return 0;
}
int D4_wavelet_2D_backward( double ** dst, double * sca_sca, double * sca_wav, double * wav_sca, double * wav_wav, int len_row, int len_col){


	
	double * hor_scale = NULL;
	double * hor_wavelet = NULL;
	
	D4_wavelet_backward_rows( &hor_scale, sca_sca, sca_wav , len_row, len_col);
	D4_wavelet_backward_rows( &hor_wavelet, wav_sca, wav_wav , len_row, len_col);
	D4_wavelet_backward_columns( dst, hor_scale, hor_wavelet, len_row, 2*len_col);

	return 0;
}

int D4_wavelet_2D_backward_CR( double ** dst, double * sca_sca, double * sca_wav, double * wav_sca, double * wav_wav, int len_row, int len_col){


	
	double * ver_scale = NULL;
	double * ver_wavelet = NULL;

	D4_wavelet_backward_columns( &ver_scale, sca_sca, sca_wav , len_row, len_col);
	D4_wavelet_backward_columns( &ver_wavelet, wav_sca, wav_wav , len_row, len_col);
	D4_wavelet_backward_rows( dst, ver_scale, ver_wavelet, 2*len_row, len_col);

	return 0;
}

int compute_dual_filter(double ** filter_dst, double * filter_src, int filter_length){
	if(filter_src == NULL || filter_dst == NULL){
		fprintf(stderr, "ERR : wavelet.c : compute_backward_filter() : At least one of the filters is a NULL pointer.\n");
		return -1;
	}
	if(filter_length < 1){
		fprintf(stdout, "ERR : wavelet.c : compute_backward_filter() : filter_length has an incorrect value : %d.\n", filter_length);
		return -1;
	}
	if(*filter_dst == NULL){
		*filter_dst = malloc(filter_length * sizeof(double));
		if(*filter_dst == NULL){
			fprintf(stderr, "ERR : wavelet.c : compute_backward_filter() : Memory allocation of filter_dst failed.\n");
			return -1;
		}
	}

	int i = 0;
	for(i = 0; i < filter_length; i++){
		(*filter_dst)[i] = filter_src[filter_length - 1 - i];
	}
	if(filter_length % 2){
		filter_length--;
	}
	double temp_array[filter_length/2];
	for(i = 0; i < filter_length / 2; i++){
		temp_array[i] = (*filter_dst)[2 * i + 1];
		(*filter_dst)[i] = (*filter_dst)[2*i];
	}
	for(i = 0; i < filter_length/2; i++){
		(*filter_dst)[ filter_length /2 + i] = temp_array[i];
	}
	return 0;
	
}


