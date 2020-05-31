#include <stdlib.h>
#include <stdio.h>

#include <cyan/image/transforms.h>

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

