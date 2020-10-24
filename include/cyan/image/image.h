#ifndef CYAN_IMAGE_H
#define CYAN_IMAGE_H

#include <cyan/common/config.h>
#include <cyan/color/color.h>

typedef struct {
	int rows ; 
	int cols ; 
    int monochrome ;                    // if set to zero, only the Y plane will be used
	double* X ;                         // XYZ as defined in CIE 1976
    double* Y ;
    double* Z ;
    enum cyan_refwhite illuminant ;     // Default to CYAN_D50
			
	size_t pixel_data_size ;
	void*  pixel_data ;
} image_t ;

image_t* image_new( int cols, int rows, int monochrome ) ;
void     image_free( image_t* ) ;
image_t * image_new_empty( int cols, int rows);

int image_cat_hor( image_t ** dst, image_t * img_left, image_t * img_right );
int image_cat_ver( image_t ** dst, image_t * img_up, image_t * img_bot );

int image_crop_rows(image_t ** dst, image_t * src, int first_row, int last_row);
int image_set_color_to_coordinates(image_t * img, color_t color, int * y_coord, int step);

int image_allocate_data_default  ( image_t*, size_t size, void* default_data ) ;
int image_allocate_data_fct( image_t*, 
                 size_t datasize, 
			     int (*fill_fct) (image_t*, int, int, void*),
			     void* context ) ;

int image_import_data   ( image_t*, size_t data_size, void* data ) ;

int image_get_data_pointer  ( image_t*, int i, int j, void** ) ;

int image_clone ( image_t*, image_t** ) ;
int image_resize ( image_t*, int cols, int rows, int monochrome, void* default_pixel_data ) ;  

int image_save(image_t* img, char* filename ) ;

image_t* image_load(char* filename , int* result ) ;

int image_strip_data(image_t*) ;


#endif
