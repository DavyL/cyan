#include <stdio.h>
#include "cyan/cyan.h"


int main( int argc, char** argv, char* envv ) {

    int ret ;
    image_t* image = NULL ; 
    
    if ( argc != 2 ) {
        fprintf( stderr, "Usage: %s <cyan image file>\n", argv[0] ) ;
        return 0 ;
    }

    ret = image_load( &image, argv[1] ) ;

    if (ret != ERR_OK ) {
        fprintf(stderr,"Could not load file %s\n", argv[1] ) ; 
        return 0 ;
    }

    printf("\tRows: %d \n", image->rows ) ;
    printf("\tCols: %d \n", image->cols ) ;
    printf("\tMonochrome: %s \n", (image->monochrome)?"yes":"no" ) ;
    printf("\tPixel Data size : %d \n", (int) image->pixel_data_size ) ;
   
    return 0 ;
}

