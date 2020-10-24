#include <string.h>
#include "cyan/common/error.h"
#include "cyan/common/config.h"
#include "cyan_file.h"


int write_header( FILE* fd, const char header[] ) {
    int ret ;
    unsigned char version[2] ;
    
    ret=fwrite( header, strlen(header), 1, fd ) ;
    if ( ret != 1 ) {
        CYAN_ERROR( ERR_FILE ) ;
        return ERR_FILE ;
    }
    
    version[0] = CYAN_VERSION_MAJOR ;
    version[1] = CYAN_VERSION_MINOR ;

    ret=fwrite( version, sizeof(unsigned char), 2, fd ) ;
    if ( ret != 2 ) {
        CYAN_ERROR( ERR_FILE ) ;
        return ERR_FILE ;
    }
    
    return ERR_OK ;
}

int read_header( FILE* fd, const char header[] ) {
    int ret ;
    char buff[10]={'\0'} ;
    unsigned char version[2] ;
    
    ret = fread(buff, strlen(header),1,fd ) ;
    if ( ret != 1 ) {
        CYAN_ERROR( ERR_FILE ) ;
        return ERR_FILE ;
    }

    if ( strcmp( buff, header ) != 0 ) {
        CYAN_ERROR( ERR_FILE ) ;
        return ERR_FILE ;
    }

    ret = fread(version, sizeof(unsigned char), 2, fd ) ;
    if ( ret != 2 ) {
        CYAN_ERROR( ERR_FILE ) ;
        return ERR_FILE ;
    }
    
    if ((version[0]!=CYAN_VERSION_MAJOR)||(version[1]!=CYAN_VERSION_MINOR)) {
        fprintf(stderr,"Cannot read file : Cyan version is not the same \n") ;
        CYAN_ERROR( ERR_FILE ) ;
        return ERR_FILE ;
    }



    return ERR_OK ;

}
