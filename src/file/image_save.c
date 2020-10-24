#include "cyan/image/image.h"
#include "cyan/common/error.h"
#include "cyan_file.h"

int image_save( image_t* img, char* filename ) {
    int ret ;
    FILE* fd = NULL ;
    
    // Open file

    fd = fopen( filename, "w" ) ;
    if ( fd == NULL ) {
        CYAN_ERROR( ERR_FILE ) ;
        return ERR_FILE ;
   }

    // write cyan header

    ret = write_header( fd, CYAN_HEADER ) ;
    if (ret){
        return ret ;
    } 

    
    // write image header

    ret = write_header( fd, CYAN_HEADER_IMAGE ) ;
    if (ret){
        return ret ;
    } 

    // write image data
    
    ret = 0 ;
    ret += fwrite( &(img->rows), sizeof(int), 1, fd ) ;
    ret += fwrite( &(img->cols), sizeof(int), 1, fd ) ;
    ret += fwrite( &(img->monochrome), sizeof(int), 1, fd ) ;
    ret += fwrite( &(img->illuminant), sizeof(enum cyan_refwhite), 1, fd ) ;
    ret += fwrite( &(img->pixel_data_size), sizeof(size_t), 1, fd ) ;

    if ( ret != 5 ) {
        CYAN_ERROR( ERR_FILE ) ;
        return ERR_FILE ;
    }

    ret = fwrite( img->Y, sizeof(double), img->rows*img->cols, fd ) ;
    if ( ret != img->rows*img->cols ) {
        CYAN_ERROR( ERR_FILE ) ;
        return ERR_FILE ;
    }

    if ( !img->monochrome ) {
        ret = fwrite( img->X, sizeof(double), img->rows*img->cols, fd ) ;
        if ( ret != img->rows*img->cols ) {
            CYAN_ERROR( ERR_FILE ) ;
            return ERR_FILE ;
        }   
        ret = fwrite( img->Z, sizeof(double), img->rows*img->cols, fd ) ;
        if ( ret != img->rows*img->cols ) {
            CYAN_ERROR( ERR_FILE ) ;
            return ERR_FILE ;
        }   
    }

    if ( img->pixel_data_size != 0 ) {
        ret = fwrite( img->pixel_data, img->pixel_data_size, img->rows*img->cols, fd ) ;
        if ( ret != img->rows*img->cols ) {
            CYAN_ERROR( ERR_FILE ) ;
            return ERR_FILE ;
        }   
    }

   // Close file

   ret=fclose(fd) ;
   if ( ret != 0 ) {
        CYAN_ERROR( ERR_FILE ) ;
        return ERR_FILE ;
   }

   return ERR_OK ;
}
