#ifndef CYAN_FILE_H
#define CYAN_FILE_H

#include <stdio.h>

#define CYAN_HEADER         "CYAN"
#define CYAN_HEADER_IMAGE   "IMG"
#define CYAN_HEADER_CLOUD   "CLD"


int write_header( FILE* fd, const char header[] ) ;
int read_header( FILE* fd , const char header[] ) ;

#endif
