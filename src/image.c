#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <cyan/image/image.h>
#include <cyan/color/color.h>

image_t* image_new(int cols, int rows) {
	image_t *tmp;
	tmp = (image_t *) malloc(sizeof(image_t));
	if (tmp == NULL) {
		fprintf(stderr, "image_new : malloc error\n");
		return NULL;
	}
	tmp->rows = rows;
	tmp->cols = cols;
    tmp->monochrome = 0 ;
    tmp->illuminant = CYAN_D50 ;
	tmp->X = (double *) malloc(rows * cols * sizeof(double));
	tmp->Y = (double *) malloc(rows * cols * sizeof(double));
	tmp->Z = (double *) malloc(rows * cols * sizeof(double));
	if ( (tmp->X == NULL) || (tmp->Y == NULL) || (tmp->Z == NULL) ) {
		fprintf(stderr, "image_new : pixels allocation error \n");
		free(tmp);
		return NULL;
	}
	tmp->pixel_data_size = (size_t) 0;
	tmp->pixel_data = (void *) NULL;
	return tmp;
}
//image_new_empty is like image_new if it had a calloc memory allocation call
image_t * image_new_empty( int cols, int rows){
	image_t * ret =NULL;
	ret = image_new( cols, rows );
	if(!ret){
		fprintf(stderr, "image_new_empty : allocation error\n");
		return NULL;
	}
	memset(ret->X, 0, rows * cols * sizeof(double));
	memset(ret->Y, 0, rows * cols * sizeof(double));
	memset(ret->Z, 0, rows * cols * sizeof(double));
	return ret;
}

void image_free(image_t * img) {
	if (img->pixel_data != (void *) NULL) {
		free(img->pixel_data);
		img->pixel_data_size = 0;
	}
	if (img->X != (double *) NULL) {
		free(img->X);
	}
	if (img->Y != (double *) NULL) {
		free(img->Y);
	}
	if (img->Z != (double *) NULL) {
		free(img->Z);
	}
	free(img);
}

//WARNING + TODO : memory mapping of image->X, Y, Z, is done the wrong way
//Concatenates two image horizontally
//img_left and img_right do not need to have the same number of rows
int image_cat_hor( image_t ** dst, image_t * img_left, image_t * img_right ){
	if(dst == NULL || img_left == NULL ||img_right == NULL){
		fprintf(stderr, "ERR : image_cat_hor : an argument is NULL.\n");
		return -1;
	}
	int rows_diff = img_left->rows - img_right->rows;
	image_t * backup = NULL;
	if(rows_diff){
		image_t * fill_image = NULL;
			
		image_t * filled_image = NULL;
		
		if(rows_diff > 0){
			fill_image = image_new_empty(img_right->cols, abs(rows_diff));
			image_cat_ver( &filled_image, img_right, fill_image);
			backup = img_right;
			img_right = filled_image;
		}else{
			fill_image = image_new_empty(img_left->cols, abs(rows_diff));
			image_cat_ver( &filled_image, img_left, fill_image);
			backup = img_left;
			img_left = filled_image;
		}
		image_free( fill_image );
	}
	
	if(*dst == NULL){
		fprintf(stdout, "image_cat_hor : dst is pointing to a NULL, allocating new image (size : %d).\n", (img_left->cols + img_right->cols) * img_left->rows * sizeof(double) );
		*dst = image_new( img_left->cols + img_right->cols, img_left->rows);
		if( *dst == NULL){
			fprintf(stderr, "ERR: image_cat_hor: Image allocation failed.\n");
			return -1;
		}	
	}


	int i = 0;

	for(i = 0; i < (*dst)->rows; i++){
		//Copying left image
		memcpy( &((*dst)->X[i*(*dst)->cols]), &(img_left->X[i * img_left->cols]), img_left->cols*sizeof(double) );
		memcpy( &((*dst)->Y[i*(*dst)->cols]), &(img_left->Y[i * img_left->cols]), img_left->cols*sizeof(double) );
		memcpy( &((*dst)->Z[i*(*dst)->cols]), &(img_left->Z[i * img_left->cols]), img_left->cols*sizeof(double) );

		memcpy( &((*dst)->X[i*(*dst)->cols + img_left->cols]), &(img_right->X[i*img_right->cols]), img_right->cols*sizeof(double)); 
		memcpy( &((*dst)->Y[i*(*dst)->cols + img_left->cols]), &(img_right->Y[i*img_right->cols]), img_right->cols*sizeof(double)); 
		memcpy( &((*dst)->Z[i*(*dst)->cols + img_left->cols]), &(img_right->Z[i*img_right->cols]), img_right->cols*sizeof(double)); 
	}
	if(rows_diff){
		if(rows_diff > 0){
			image_free(img_right);
			img_right = backup;
		}else{
			image_free(img_left);
			img_left = backup;
		}
	}


	fprintf(stdout, "image_cat_hor : Success\n");

	return 0;

}

int image_cat_ver( image_t ** dst, image_t * img_up, image_t * img_bot ){
	if(dst == NULL || img_up == NULL ||img_bot == NULL){
		fprintf(stderr, "ERR : image_cat_ver : an argument is NULL.\n");
		return -1;
	}
	int cols_diff = img_up->cols - img_bot->cols;
	image_t * backup = NULL;
	if(cols_diff){
		image_t * fill_image = NULL;	
		image_t * filled_image = NULL;
		
		if(cols_diff > 0){
			fill_image = image_new_empty(abs(cols_diff), img_bot->rows);
			image_cat_hor( &filled_image, img_bot, fill_image);
			backup = img_bot;
			img_bot = filled_image;
		}else{
			fill_image = image_new_empty(abs(cols_diff), img_up->rows);
			image_cat_hor( &filled_image, img_up, fill_image);
			backup = img_up;
			img_up = filled_image;
		}
		image_free( fill_image );
	}
	
	
	
	if(*dst == NULL){
		fprintf(stdout, "image_cat_ver : dst is pointing to a NULL, allocating new image (size : %d).\n", (img_up->cols) * (img_up->rows + img_bot->rows) * sizeof(double) );
		*dst = image_new( img_up->cols, img_up->rows + img_bot->rows);
		if( *dst == NULL){
			fprintf(stderr, "ERR: image_cat_ver: Image allocation failed.\n");
			return -1;
		}	
	}

	int i = 0;
	for(i = 0; i < img_up->rows; i++){

		//Copying upper image
		memcpy( &((*dst)->X[i*(*dst)->cols]), &(img_up->X[i * img_up->cols]), img_up->cols*sizeof(double)  );
		memcpy( &((*dst)->Y[i*(*dst)->cols]), &(img_up->Y[i * img_up->cols]), img_up->cols*sizeof(double)  );
		memcpy( &((*dst)->Z[i*(*dst)->cols]), &(img_up->Z[i * img_up->cols]), img_up->cols*sizeof(double)  );
	}
	for( i = 0; i < img_bot->rows; i++ ){

		//Copying bottom image
		memcpy( &((*dst)->X[(i + img_up->rows)*(*dst)->cols]), &(img_bot->X[i*img_bot->cols]), img_bot->cols*sizeof(double)); 
		memcpy( &((*dst)->Y[(i + img_up->rows)*(*dst)->cols]), &(img_bot->Y[i*img_bot->cols]), img_bot->cols*sizeof(double)); 
		memcpy( &((*dst)->Z[(i + img_up->rows)*(*dst)->cols]), &(img_bot->Z[i*img_bot->cols]), img_bot->cols*sizeof(double)); 
	}
	if(cols_diff){
		if(cols_diff > 0){
			image_free(img_bot);
			img_bot = backup;
		}else{
			image_free(img_up);
			img_up = backup;
		}
	}


	fprintf(stdout, "image_cat_ver : Success\n");

	return 0;	
}

//"Crops" an image horizontally
//Stores in dst an image with last_row - first_row rows and src->cols columns
//image_crop_rows doesn't make a copy, data is the same in dst and src, except that the data pointed by dst has an offset from src
//with the appropriate change of value for (*dst)->rows
//It's only a zoomed image
//Do not use free_image on dst, use it on src instead (just use free() on dst)
int image_crop_rows(image_t ** dst, image_t * src, int first_row, int last_row){
	
	int new_rows = last_row - first_row;
	if( dst == NULL || src == NULL || new_rows <= 0 || first_row < 0){
		fprintf(stdout, "image_crop_rows : NULL pointer or unappropriate row values were given as arguments.\n");
		return -1;
	}
	if( last_row > src->rows){
		fprintf(stdout , "image_crop_rows : number of rows to crop is larger than the number of rows in src.\n");
		return -1;
	}
	*dst = (image_t * )malloc(sizeof(image_t) );
	if(*dst == NULL){
		fprintf(stderr, "image_crop_row : allocation error\n");
		return -1;
	}
	(*dst)->cols = src->cols;
	(*dst)->rows = new_rows;
	
	(*dst)->pixel_data_size = src->pixel_data_size;
	//An error might happen here if pixel_data_size is not zero and not enough (or no space at all) has been allocated to data
	//Hence, we check if src->pixel_data is NULL to prevent errors 
	//If pixel_data has not been initialized to NULL, the function is not safe
	if( !src->pixel_data){
		(*dst)->pixel_data = src->pixel_data + first_row * src->cols * (*dst)->pixel_data_size;
	}else{
		(*dst)->pixel_data = NULL;
	}
	(*dst)->monochrome = src->monochrome;

	(*dst)->X = src->X + first_row * src->cols;
	(*dst)->Y = src->Y + first_row * src->cols;
	(*dst)->Z = src->Z + first_row * src->cols;
	
	
	return 0;
}
int image_set_color_to_coordinates(image_t * img, color_t color, int * y_coord, int step){
	
	if(img == NULL || y_coord == NULL){
		fprintf(stderr, "image_set_color_to_coordinates: a NULL pointer was passed\n");
		return -1;
	}
	if(step <= 0){
		fprintf(stderr, "image_set_color_to_coordinates: step doesn't have a correct value (<= 0)\n");
		return -1;
	}

	int i = 0;
	int temp_y_coord;
	for( i = 0; i < img->cols - step; i+= step){
		temp_y_coord = y_coord[i];
		//following conditions to ensure values are in range
		if(temp_y_coord > img->rows){
			temp_y_coord = img->rows - 1;
		}else if(temp_y_coord < 0){
			temp_y_coord = 0;
		}
		img->X[ img->cols*(img->rows - temp_y_coord ) + i ] = color.X;
		img->Y[ img->cols*(img->rows - temp_y_coord ) + i ] = color.Y;
		img->Z[ img->cols*(img->rows - temp_y_coord ) + i ] = color.Z;
	}
	return i;
}

int image_allocate_data_default(image_t * img, size_t size, void *data_array) {
	int i;
	if (img == (image_t *) NULL) {
		fprintf(stderr,
			"image_allocate_data_default: image is a NULL pointer\n");
		return -1;
	}
	if (img->pixel_data != (void *) NULL) {
		fprintf(stderr,
			"image_allocate_data_default : img->pixel_data is not NULL \n");
		return -1;
	}
	img->pixel_data_size = size;
	img->pixel_data = (void *) malloc(img->rows * img->cols * size);
	if (img->pixel_data == NULL) {
		fprintf(stderr,
			"image_allocate_data_default : data allocation error \n");
		img->pixel_data_size = 0;
		return 0;
	}
	for (i = 0; i < img->rows * img->cols; i++) {
		memcpy(img->pixel_data + i * size, data_array, size);
	}
	return i;
}

int image_allocate_data_fct(image_t * img, size_t size, int (*fill_fct) (image_t *, int, int, void *), void *context) {
	int i, j, ret;
	if (img == (image_t *) NULL) {
		fprintf(stderr,
			"image_allocate_data_fct: image is a NULL pointer\n");
		return -1;
	}
	if (img->pixel_data != (void *) NULL) {
		fprintf(stderr,
			"image_allocate_data_fct : img->pixel_data is not NULL \n");
		return -1;
	}
	img->pixel_data_size = size;
	img->pixel_data = (void *) malloc(img->rows * img->cols * size);
	if (img->pixel_data == NULL) {
		fprintf(stderr,
			"image_allocate_data_fct : data allocation error \n");
		img->pixel_data_size = 0;
		return 0;
	}
	ret = 0;
	for (j = 0; j < img->rows; j++)
		for (i = 0; i < img->cols; i++) {
			ret += fill_fct(img, i, j, context);
		}
	return ret;
}

int image_import_data(image_t * img, size_t data_size, void *data) {
	if (img == (image_t *) NULL) {
		fprintf(stderr,
			"image_import_data: image is a NULL pointer\n");
		return -1;
	}
	if (img->pixel_data != (void *) NULL) {
		fprintf(stderr,
			"image_import_data : img->pixel_data is not NULL \n");
		return -1;
	}
	img->pixel_data_size = data_size;
	memcpy(img->pixel_data, data, img->cols * img->rows * data_size);
	return img->cols * img->rows;
}

void *image_get_data_pointer(image_t * img, int i, int j) {
	return img->pixel_data + (i + j * img->cols) * img->pixel_data_size;
}

image_t *image_clone(image_t * img) {
	image_t *clone;
	clone = image_new(img->cols, img->rows);
    clone->monochrome = img->monochrome ;
    clone->illuminant = img->illuminant ;
    fprintf(stdout, "Inside image_clone().\n");
    memcpy(clone->Y, img->Y, img->rows * img->cols * sizeof(double));
	if ( img->monochrome ) {
		//HOTFIX
		//There should be a better way to do this (i.e, without allocating a
		//whole array of zeros which should not be used 
        	//Assigning a NULL pointer to either of these causes segmentation fault upon allocation after exiting image_clone()
	clone->X = calloc(img->rows*img->cols, sizeof(double)) ;
        clone->Z = clone->X;
    } else {
	    memcpy(clone->X, img->X, img->rows * img->cols * sizeof(double));
	    memcpy(clone->Z, img->Z, img->rows * img->cols * sizeof(double));
    }
    if (img->pixel_data != (void *) NULL) {
		image_import_data(clone, img->pixel_data_size, img->pixel_data);
	}
    fprintf(stdout, "About to exit image_clone().\n");
	return clone;
}

