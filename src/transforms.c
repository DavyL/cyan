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

//FT_image_Y allows to compute FFT_2D of image in the direction specified by transform (FFT_2D or FFT_2D_reverse should be given as an argument)
//The coefficients computed and saved in ft_array are given as parameters to convert_cart_to_image
//which specifies how coefficients should be stored in the image
//FT_image_Y may actually be used by any kind of transforms which respects inputs and output types
image_t * FT_image_Y(image_t * image, complex_cart_t ** (*transform)(complex_cart_t **, int, int), image_t * (*convert_cart_to_image)( complex_cart_t **, int, int) ){
 	if(image == NULL)
		fprintf(stderr, "FT_image : image is a NULL pointer. \n");

	//1. Convert image in array[rows][cols]
	complex_cart_t ** image_array = NULL;
	image_array = image_to_cart(image, Y_to_cart); 
	//2. Compute 2D FFT
	complex_cart_t ** ft_array = NULL;
	ft_array = transform(image_array, 9, 9);
	//3. Convert 2D FFT in an image
	image_t * image_ft = NULL;
	
	image_ft = convert_cart_to_image(ft_array, image->cols, image->rows);

	return image_ft;
}
image_t * cart_to_Y(complex_cart_t ** array, int N, int M){
	return cart_to_image(NULL, array, N, M, zero, sum, zero);
}

image_t * cart_power_to_Y(complex_cart_t ** array, int N, int M){
	return cart_to_image(NULL, array, N, M, zero, power_complex_cart, zero);
}
image_t * cart_log_power_to_Y(complex_cart_t ** array, int N, int M){
	return cart_to_image(NULL, array, N, M, zero, log_power_complex_cart, zero);
}
image_t * cart_phase_to_Y(complex_cart_t ** array, int N, int M){
	return cart_to_image(NULL, array, N, M, zero, phase_complex_cart, zero);
}
 

complex_cart_t ** image_to_cart(image_t * image, complex_cart_t (*xyz_to_cart)( double, double, double) ){
	if(image == NULL)
		fprintf(stderr, "image_to_cart : image is a NULL pointer. \n");

	complex_cart_t ** array_cart = malloc( image->rows * sizeof(complex_cart_t *));;

	if(array_cart == NULL){
		fprintf(stderr, "image_to_cart : couldn't allocate memory to array_cart. \n");
	}
	int i, j, coords;
	for(i = 0; i < image->rows; i++){
		array_cart[i] = malloc(image->cols * sizeof(complex_cart_t));
		for(j=0; j<image->cols; j++){
			coords = j + i * image->cols;
			array_cart[i][j] = xyz_to_cart(image->X[coords], image->Y[coords], image->Z[coords] );
		}
	}
	return array_cart;
}
image_t * cart_to_image(image_t * image, complex_cart_t ** array, int N, int M,
	       double (*X_func)(complex_cart_t), double (*Y_func)(complex_cart_t), double (*Z_func)(complex_cart_t)	){
	if(array == (complex_cart_t **) NULL)
		fprintf(stderr, "cart_to_image : array is a NULL pointer. \n");

	image_t * ret;
	if(image == (image_t *) NULL)
		ret = image_new(N, M);
	else
		ret = image;

	//complex_cart_t ** array_cart = NULL;
	//fprintf(stdout, "rows : %d, cols : %d\n", image->rows, image->cols);
	int i, j, coord;
	for(i = 0; i < N; i++){
		for(j=0; j < M; j++){
			coord = j + i * M;
			ret->X[coord] = X_func(array[i][j]) ;
			ret->Y[coord] = Y_func(array[i][j]) ;
			ret->Z[coord] = Z_func(array[i][j]) ;
		}
	}
	return ret;
}
complex_cart_t X_Y_Z_to_cart(double X, double Y, double Z){
	complex_cart_t cart;

	cart.real =X + Y + Z ;
	//cart.im = X + Z;

	return cart;

}
complex_cart_t X_Y_to_cart(double X, double Y, double Z){
	complex_cart_t cart;

	cart.real = Y;
	cart.im = X;

	return cart;

}
complex_cart_t  Y_to_cart(double X, double Y, double Z){
	complex_cart_t  cart;
	//cart = (complex_cart_t *) malloc(sizeof(complex_cart_t));
	cart.real = Y;
	cart.im =0.0f;	

	return cart;
}
//lowpass and highpass are temporary functions
//They shall be properly developped later
complex_cart_t ** lowpass_harsh(complex_cart_t ** array_cart, int N, int M){
	complex_cart_t ** ret = NULL;
	ret = malloc( N * sizeof(complex_cart_t *) );
	int i,j;
	for(i = 0; i < N; i++){
		ret[i] = malloc(M*sizeof(complex_cart_t));

		if( i < N/4 ||i > 3*N/4){
			for(j = 0; j < M; j++){		
				ret[i][j].real = 0.0f;
				ret[i][j].im = 0.0f;
			}
		}else{
			for(j = 0; j < M; j++){		
				ret[i][j].real = array_cart[i][j].real;
				ret[i][j].im = array_cart[i][j].im;
			}
		}
	}
	return ret;
}
complex_cart_t ** middle_band_pass( complex_cart_t ** array_cart, int N, int M){

	complex_cart_t ** ret = NULL;
	ret =(complex_cart_t **)  malloc( N * sizeof(complex_cart_t *));

	int i, j;
	for(i = 0; i < N; i++){
		ret[i] = ( complex_cart_t * ) malloc(M * sizeof(complex_cart_t) );
			for(j = 0; j< M; j++){
				if( !(j > M /4 && j < 3*M/4)){ 
					ret[i][j].real = array_cart[i][j].real;
					ret[i][j].im = array_cart[i][j].im;
				}else{
					ret[i][j].real = 0.0f;
					ret[i][j].im = 0.0f;
				}
			}

	}
	return ret;
}

int rotate_180( void *** dst, void ** src, int N, size_t size){
	int i, j;
	if(*dst ==NULL){
		*dst = (void ** ) malloc( N * sizeof(void *));
		if(*dst == NULL){
			fprintf( stderr, "rotate_180 : memory allocation failed. \n" );
			return -1;
		}
	}
	for( i = 0; i < N; i++){
		memcpy( *dst + i, src + N - 1 - i, size );
	}

	return 0;
}
int transpose( void *** dst, void ** src, int N, size_t size ){
	//src is a N elem array, each elem being a pointer to an N elem array of elem of size "size"
	//Transpose sets each line to be a column
	//dst will be an N elem array, each elem being a pointer to an N elem array of elem of size "size"
	//dst is pointing to the array, hence dst shall be given as the address of a pointer in order to allow memory allocation if pointing to NULL
	//If dst is not a NULL pointer and pointing to a memory area which is not of the proper size or type, errors or unexpected behaviour could happen 

	//Note that transpose only works for squared arrays, a version which fills with zero the undefined elements could be defined in the future	
	int i, j;
	if(*dst == NULL){
		*dst = (void **)  malloc( N * sizeof( void *) );
	} 
	for(i = 0; i < N; i++){
		(*dst)[i] = malloc( N * size );
		if( (*dst)[i] == NULL){
			fprintf(stderr, "transpose : memory allocation failed.\n");
			return -1;

		}
		for(j = 0; j < N; j++){
			memcpy( *(*dst + i) + j * size,  *( src + j) + i * size, size );
		}
	}
	return 0;
}

void * rotate_buffer( void * buffer, int N, size_t size){
	
	void * temp = malloc(N/2 * size);

	memcpy(temp, buffer + (N/2 - 1) * size, N/2 * size);	//temp contains the second half of buffer
	
	memcpy(buffer + (N/2 -1) * size, buffer, N/2 * size);	//The second half of buffer now contains the previous first half
	memcpy(buffer, temp, N/2 * size);	//And its first half is filled with temp

	free(temp);

	return buffer;
}
//reverse buffer maps buffer[k] ---> buffer[N-k]
void * reverse_buffer( void * buffer, int N, size_t size){
	void * temp = malloc(N/2 * size);
	memcpy(temp, buffer, N/2 * size);

	int i;
	for( i=0; i< N/2; i++){
		memcpy(buffer + i * size, buffer + (N - i - 1)*size, size);	//Shouldn't we have a -1 in here ? 
									//I would suspect to have an out of bound here 
									//but it seems to work this way...
		memcpy( buffer + (N - i - 1) * size, temp + i * size, size);
	}

	free(temp);

	return buffer;
}

	
//Returns a grey level copy of color_img
image_t * color2grey(image_t * color_img){

	fprintf(stdout, "Entering color2grey().\n");

	if(color_img == NULL){
		fprintf(stderr, "color2grey() was given a NULL pointer instead of an image.\n");
		return NULL;
	}
	
	//Saving monochrome state of color_img	
	int monochrome_backup = color_img->monochrome;
	color_img->monochrome = 1;
	
	//Making a copy of color_img w/ a monochrome state
	
	fprintf(stdout, "Cloning color_img.\n");
//	fprintf(stdout, "clone : %x, image : %x\n",image_clone( color_img ), color_img);
	
	image_t * grey_img;
	grey_img = image_clone(color_img);	

	//Restoring monochrome state in color_img in order to keep it unchanged through color2grey()
	color_img->monochrome = monochrome_backup;

	fprintf(stdout, "About to exit color2grey()");
	return grey_img;

}

//dot_product_4_f() returns the dot product of 2 vectors containing 4 floats
float dot_product_4_f(	float x1, float x2, float x3, float x4,
			float y1, float y2, float y3, float y4 ){
	return x1*y1 + x2*y2 + x3*y3 + x4*y4 ;
}

//Returns the norm 2 (euclidean) squared of a vector containing 4 floats
float norm_2_sq_4_f( float x, float y, float z, float t){
	return dot_product_4_f( x, y, z, t, x, y, z, t) ;
}
//Sets to 0 every value which is absolutely smaller than T
//Returns -1 if there is an error, otherwise it returns the number of values kept in the array
int double_threshold( double * array, double T, int len){
	if(array == NULL){
		fprintf(stderr, "Invalid argument : double_threshold() arrray is a NULL pointer.\n");
		return -1;
	}
	int i = 0;
	for( i = len - 1; i >= 0; i--){
		if(fabs(array[i]) < T){
			array[i] = 0.0f;
			len--;
		}
	}
	return len;
}
//Compute the l2 distance of two arrays
int l2_distance(double * dst, double * x, double * y, int len){
	if(dst == NULL || x == NULL){
		fprintf(stderr, "Invalid argument: l2_distance() a NULL pointer was passed as an argument.\n");
		return -1;
	}
	
	*dst = 0.0f;
	int i = 0;	
	if( y != NULL){
		for(i = 0; i < len; i++){
			*dst += fabs(x[i]*x[i] - y[i]*y[i]);
		}
	}
	else{
		for(i = 0; i < len; i++){
			*dst += x[i]*x[i];
		}
	}
	*dst = sqrt(*dst);
	return 0;


}

//Stores in a the convolution of 1D array (of length n) (which is a real value, i.e. a double)
//at step step with a filter h of length filter_len
//Warning:  a is not initialized
int convol_loc_1D( double * a, double * array, double * h, int step, int filter_len, int len){
	if( a == NULL || array == NULL || h == NULL){
		fprintf(stderr, "ERR : transforms.c : convol_loc_1D() : a NULL pointer was passed as an argument.\n");
		return -1;
	}
	if(step < - filter_len/2 || step >= len){
		fprintf(stderr, "ERR : transforms.c : convol_loc_1D() : an out of bound value is passed as an argument for step : %d.\n", step);
		return -1;
	}

	if(filter_len % 2 != 0){		//TODO : Implement this later
		fprintf(stderr, "ERR : transforms.c : convol_loc_1D() : filter_len is an odd value.\n");
		return -1;
	}
		
	//Ensuring the first element is in range
	int first_elem = 0;	
	if( step - filter_len /2 >= 0){
		first_elem = step - filter_len / 2;
	}
	
	
	int j = 0;
	for( j = 0; j < filter_len; j++){
		if(first_elem + j < 0){
			*a += array[ - first_elem - j] * h[j];
		}else if( first_elem + j < len){
			//*a =*a + (*(array + first_elem + j )) * (h[j]);
			*a += array[ first_elem + j] * (h[j]);
		}else{

			*a =*a + (*(array + first_elem + filter_len - j )) * (h[j]);	
			//Making symmetry on the boundary
			//Modify here to change boundaries behaviour
		}
	}
	
	return 0;
}
//Performs the convolution on columns from array and stores the result as a double in a
//The input array is of size len_row*len_col
//Hence the offset between each element on which to perform the convolution is len_col
//Location input is (row, column) the position of the first elem
int convol_loc_1D_columns( double * a, double * array, double * h, int row, int column, int filter_len, int len_row, int len_col){
	if( a == NULL || array == NULL || h == NULL){
		fprintf(stderr, "ERR : transforms.c : convol_loc_1D_columns() : a NULL pointer was passed as an argument.\n");
		return -1;
	}
	if(row < 0 || row >= len_row){
		fprintf(stderr, "ERR : transforms.c : convol_loc_1D_columns() : an out of bound value is passed as an argument for i : %d.\n", row);
		return -1;
	}

	if(filter_len % 2 != 0){		//TODO : Implement this later
		fprintf(stderr, "ERR : transforms.c : convol_loc_1D() : filter_len is an odd value.\n");
		return -1;
	}
		
	//Ensuring the first element is in range
	int first_elem = 0;	//TODO : Add "wrapping" on the first boundary in the 'for' loop 
	if( row - filter_len /2 >= 0){
		first_elem = row - filter_len / 2;
	}
	
	
	int j = 0;
	for( j = 0; j < filter_len; j++){
		if( first_elem + j < len_row){
			*a += array[ (first_elem + j)*len_col + column] * (h[j]);
		}
		else{

			*a =*a + array[ (first_elem + filter_len - j)*len_col + column] * (h[j]);	
			//Making symmetry on the boundary
			//Modify here to change boundaries behaviour
		}
	}
	
	return 0;
}

//Performs full 2D convolution of an array
int convol_2D( double ** dst, double * array, int array_rows, int array_cols, double * matrix, int matrix_rows, int matrix_cols){
	if(dst == NULL || array == NULL || matrix == NULL){
		fprintf(stderr, "ERR : transforms.c : convol_2D() : One of the arguments is a NULL pointer.\n");
		return -1;
	}
	if( array_rows <= 0 || array_cols <= 0 ){
		fprintf(stderr, "ERR : transforms.c : convol_2D() : Array rows or cols is a negative number.\n");
		return -1;
	}
	if( matrix_rows <= 0 || matrix_cols <= 0 ){
		fprintf(stderr, "ERR : transforms.c : convol_2D() : Matrix rows or cols is a negative number.\n");
		return -1;
	}

	if(*dst == NULL){
		*dst = calloc(array_rows * array_cols,  sizeof(double));
		if(*dst == NULL){
			fprintf(stderr, "ERR : transforms.c : convol_2D() : Memory allocation of dst failed.\n");
			return -1;
		}
	}

	int i = 0;
	int j = 0;
	for( i = 0; i < array_rows; i++){
		for( j = 0; j < array_cols; j++){
			if(convol_2D_loc( &((*dst)[i * array_cols + j]), array, i, j, array_rows, array_cols, matrix, matrix_rows, matrix_cols) == -1){
				fprintf(stderr, "ERR : transforms.c : convol_2D() : Something went wrong at step (i,j) : (%d, %d).\n", i, j);
				return -1;
			}
		}
	}
	
	return 0;
}

//Computes matrix convolution on a local support
int convol_2D_loc( double * dst, double * array, int first_row, int first_col, int array_rows, int array_cols, double * matrix, int matrix_rows, int matrix_cols){
	if(dst == NULL){
		fprintf(stderr, "ERR : transforms.c : convol_2D_loc() : dst is a NULL pointer. It should be the address of a double variable.\n");
		return -1;
	}
	
	if( array == NULL || matrix == NULL){
		fprintf(stderr, "ERR : transforms.c : convol_2D_loc() : One of the arguments is a NULL pointer.\n");
		return -1;
	}
	if( array_rows <= 0 || array_cols <= 0 ){
		fprintf(stderr, "ERR : transforms.c : convol_2D_loc() : Array rows or cols is a negative number.\n");
		return -1;
	}
	if( matrix_rows <= 0 || matrix_cols <= 0 ){
		fprintf(stderr, "ERR : transforms.c : convol_2D_loc() : Matrix rows or cols is a negative number.\n");
		return -1;
	}

	int i = 0;
	int first_step = first_row - matrix_rows /2;
	for( i = 0; i < matrix_rows; i++){
		if(first_step + i < 0){
			convol_loc_1D(dst, array - (first_step + i) * array_cols, matrix + i * matrix_cols, first_col, matrix_cols, array_cols);
		}else if(first_step + i >= array_rows){
			convol_loc_1D(dst, array + (array_rows - i) * array_cols, matrix + i * matrix_cols, first_col, matrix_cols, array_cols);
		}else{
			convol_loc_1D(dst, array + (first_step + i) * array_cols, matrix + i * matrix_cols, first_col, matrix_cols, array_cols);
		}

	}
	
	return 0;
}

//Stores tensor product of h and g (two vectors)
int tensor_double( double ** dst, double * h, int h_len, double * g, int g_len){
	if(dst == NULL || h == NULL){
		fprintf(stderr, "ERR : transforms.c : tensor_double() : One of the arguments is a NULL pointer.\n");
		return -1;
	}
	if(*dst == NULL){
		*dst = malloc( h_len * g_len * sizeof(double));
		if(*dst == NULL){
			fprintf(stderr, "ERR : transforms.c : tensor_double() : Memory allocation failed.\n");
		}
	}
	int i = 0;
	int j = 0;
	for(i = 0; i < h_len; i++){
		for(j=0; j < g_len; j++){
			(*dst)[ i*h_len + j ] = h[i] * g[j];
		}
	}
	return 0;

}



