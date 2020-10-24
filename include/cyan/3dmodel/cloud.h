#ifndef CYAN_CLOUD_H
#define CYAN_CLOUD_H

#include <cyan/common/config.h>
#include <cyan/color/color.h>

typedef struct cloud_s {
	
	int nb_points ;
	
	double*  coords ;       // 3D Coords in XYZ order
	double*  color_X ;      // Color in XYZ colorspace as defined in CIE 1976
	double*  color_Y ;
	double*  color_Z ;

    enum cyan_refwhite illuminant ;     // Default to CYAN_D50

	size_t  point_data_size ;
	void*   point_data  ;

} cloud_t ;

cloud_t* cloud_new( int nb_points ) ;   
void cloud_free( cloud_t* ) ; 

int cloud_allocate_color_default ( cloud_t*,
                                    enum cyan_refwhite illuminant, 
                                    double X, double Y, double Z ) ;
int cloud_allocate_color_fct( cloud_t*,                                  
                                enum cyan_refwhite illuminant,
                                int (*color_fct) (cloud_t*, int, void* ),
                                void* context ) ;

int cloud_allocate_data_default  ( cloud_t*, size_t size,  void* default_data ) ;   
int cloud_allocate_data_fct ( cloud_t*, 
                                size_t datasize,
                                int (*fill_fct) ( cloud_t*, int, void* ),
                                void* context ) ;

int cloud_import_data( cloud_t*, size_t data_size, void* data ) ; 
int cloud_get_data_pointer( cloud_t*, int i, void** ) ;            

int cloud_clone ( cloud_t*, cloud_t** ) ;
int cloud_merge ( cloud_t*, cloud_t* ) ; 


int cloud_add_point( cloud_t* cloud,                               
                        double x, double y, double z, 
                        double color_x, double color_y, double color_z,
                        void* data ) ; 

int cloud_del_point ( cloud_t* cloud, int point ) ; 

void cloud_strip_color ( cloud_t*) ; 
void cloud_strip_data   ( cloud_t*) ; 

#endif
