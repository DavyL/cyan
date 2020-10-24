#include <stdlib.h>

#include "cyan/image/image.h"
#include "cyan/common/error.h"
#include "cyan_file.h"

int image_load( image_t** img, char* filename ) {
    
    int ret ;
    FILE* fd = NULL ;
    int rows ;
    int cols ;
    int monochrome ;
    enum cyan_refwhite illuminant ;
    size_t pixel_data_size ;

    // Open file

    fd = fopen( filename, "r" ) ;
    if ( fd == NULL ) {
        CYAN_ERROR( ERR_FILE ) ;
        return ERR_FILE ;
   }

    // read cyan header

    ret = read_header( fd, CYAN_HEADER ) ;
    if (ret){
        return ret ;
    } 

    
    // read image header

    ret = read_header( fd, CYAN_HEADER_IMAGE ) ;
    if (ret){
        return ret ;
    } 

    // read image data
    
    ret = 0 ;
    ret += fread( &(rows), sizeof(int), 1, fd ) ;
    ret += fread( &(cols), sizeof(int), 1, fd ) ;
    ret += fread( &(monochrome), sizeof(int), 1, fd ) ;
    ret += fread( &(illuminant), sizeof(enum cyan_refwhite), 1, fd ) ;
    ret += fread( &(pixel_data_size), sizeof(size_t), 1, fd ) ;
    
    if ( ret != 5 ) {
        CYAN_ERROR( ERR_FILE ) ;
        return ERR_FILE ;
    }

    if ( *img == NULL ) {
        *img = image_new( cols, rows, monochrome ) ;
    } else {
        image_strip_data(*img) ;
        if (((*img)->rows!=rows)||((*img)->cols!=cols)||((*img)->monochrome!=monochrome)) {
            image_resize( *img, cols, rows, monochrome, NULL ) ;
        }
    }

    (*img)->illuminant = illuminant ;
    (*img)->pixel_data_size = pixel_data_size ;

    if ( pixel_data_size != 0 ) {
        (*img)->pixel_data = malloc( rows*cols*(*img)->pixel_data_size ) ;
    } else {
        (*img)->pixel_data = NULL ;
    }

    
    ret = fread( (*img)->Y, sizeof(double), (*img)->rows*(*img)->cols, fd ) ;
    if ( ret != (*img)->rows*(*img)->cols ) {
        CYAN_ERROR( ERR_FILE ) ;
        return ERR_FILE ;
    }

    if ( !(*img)->monochrome ) {
        ret = fread( (*img)->X, sizeof(double), (*img)->rows*(*img)->cols, fd ) ;
        if ( ret != (*img)->rows*(*img)->cols ) {
            CYAN_ERROR( ERR_FILE ) ;
            return ERR_FILE ;
        }   
        ret = fread( (*img)->Z, sizeof(double), (*img)->rows*(*img)->cols, fd ) ;
        if ( ret != (*img)->rows*(*img)->cols ) {
            CYAN_ERROR( ERR_FILE ) ;
            return ERR_FILE ;
        }   
    }

    if ( (*img)->pixel_data_size != 0 ) {
        ret = fread( (*img)->pixel_data, (*img)->pixel_data_size, (*img)->rows*(*img)->cols, fd ) ;
        if ( ret != (*img)->rows*(*img)->cols ) {
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

