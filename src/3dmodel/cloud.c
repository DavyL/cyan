#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "cyan/common/error.h"
#include "cyan/3dmodel/cloud.h"

cloud_t* cloud_new( int nb_points ) { 

    cloud_t* tmp ;

	tmp = (cloud_t *) malloc(sizeof(cloud_t));
	if (tmp == NULL) {
		CYAN_ERROR( ERR_MALLOC ) ;
		return NULL;
	}

    tmp->nb_points = nb_points ;
    tmp->coords = (double*) malloc ( nb_points*3*sizeof(double) ) ;
    if ( tmp->coords == NULL ) {
        CYAN_ERROR( ERR_MALLOC ) ;
        return NULL ;
    }

    tmp->color_X = NULL ;
    tmp->color_Y = NULL ;
    tmp->color_Z = NULL ;
    tmp->illuminant = CYAN_D50 ;
    tmp->point_data_size = 0 ;
    tmp->point_data = NULL ;

    return tmp ;

}

void cloud_free( cloud_t* cloud ) {
    
    if ( cloud->coords != NULL )
        free ( cloud->coords ) ;
    
    if ( cloud->color_X != NULL )
        free ( cloud->color_X ) ;
    
    if ( cloud->color_Y != NULL ) 
        free ( cloud->color_Y ) ;
    
    if ( cloud->color_Z != NULL )
        free ( cloud->color_Z ) ;

}

int cloud_allocate_color_default(cloud_t * cld, enum cyan_refwhite illuminant, double X, double Y, double Z) {
    int i;
	if (cld == (cloud_t *) NULL) {
        CYAN_ERROR( ERR_NULL ) ;
		return ERR_NULL ;
	}
    if ((cld->color_X != (void *) NULL) ||
        (cld->color_Y != (void *) NULL) || 
        (cld->color_Z != (void *) NULL) ){
		CYAN_ERROR( ERR_INVALID_ARG ) ;
        return ERR_INVALID_ARG ;
	}
	cld->color_X = (double *) malloc(cld->nb_points * sizeof(double));
	if (cld->color_X == NULL) {
		CYAN_ERROR( ERR_MALLOC ) ;
        return ERR_MALLOC ;
	}
	cld->color_Y = (double *) malloc(cld->nb_points * sizeof(double));
	if (cld->color_Y == NULL) {
		CYAN_ERROR( ERR_MALLOC ) ;
        return ERR_MALLOC ;
	}
	cld->color_Z = (double *) malloc(cld->nb_points * sizeof(double));
	if (cld->color_Z == NULL) {
		CYAN_ERROR( ERR_MALLOC ) ;
        return ERR_MALLOC ;
	}
    cld->illuminant = illuminant ;
    for (i=0; i<cld->nb_points;i++) {
        cld->color_X[i] = X ;
        cld->color_Y[i] = Y ;
        cld->color_Z[i] = Z ;
    }
    return ERR_OK ;
}


int cloud_allocate_color_fct(cloud_t * cld, enum cyan_refwhite illuminant,
                                int (*color_fct)(cloud_t*, int,void*),
                                void* context ) {
    int i;
    int ret ;
	if (cld == (cloud_t *) NULL) {
        CYAN_ERROR( ERR_NULL ) ;
		return ERR_NULL ;
	}
    if ((cld->color_X != (void *) NULL) ||
        (cld->color_Y != (void *) NULL) || 
        (cld->color_Z != (void *) NULL) ){
		CYAN_ERROR( ERR_INVALID_ARG ) ;
        return ERR_INVALID_ARG ;
	}
	cld->color_X = (double *) malloc(cld->nb_points * sizeof(double));
	if (cld->color_X == NULL) {
		CYAN_ERROR( ERR_MALLOC ) ;
        return ERR_MALLOC ;
	}
	cld->color_Y = (double *) malloc(cld->nb_points * sizeof(double));
	if (cld->color_Y == NULL) {
		CYAN_ERROR( ERR_MALLOC ) ;
        return ERR_MALLOC ;
	}
	cld->color_Z = (double *) malloc(cld->nb_points * sizeof(double));
	if (cld->color_Z == NULL) {
		CYAN_ERROR( ERR_MALLOC ) ;
        return ERR_MALLOC ;
	}
    cld->illuminant = illuminant ;
    ret = ERR_OK ;
    for (i=0; i<cld->nb_points;i++) {
        ret += color_fct( cld, i, context ) ;
    }
    return ret ;
}

int cloud_allocate_data_default(cloud_t * cld, size_t size, void *default_data) {
	int i;
	if (cld == (cloud_t *) NULL) {
        CYAN_ERROR( ERR_NULL ) ;
		return ERR_NULL ;
	}
	if (cld->point_data != (void *) NULL) {
		CYAN_ERROR( ERR_INVALID_ARG ) ;
        return ERR_INVALID_ARG ;
	}
	cld->point_data_size = size;
	cld->point_data = (void *) malloc(cld->nb_points * size);
	if (cld->point_data == NULL) {
		CYAN_ERROR( ERR_MALLOC ) ;
        return ERR_MALLOC ;
	}
	for (i = 0; i < cld->nb_points; i++) {
		memcpy(cld->point_data + i * size, default_data, size);
	}
	return ERR_OK ;
}


int cloud_allocate_data_fct(cloud_t * cld, size_t size, 
                            int (*fill_fct) (cloud_t*, int, void*),
                            void * context) {
	int i;
    int ret ;
	if (cld == (cloud_t *) NULL) {
        CYAN_ERROR( ERR_NULL ) ;
		return ERR_NULL ;
	}
	if (cld->point_data != (void *) NULL) {
		CYAN_ERROR( ERR_INVALID_ARG ) ;
        return ERR_INVALID_ARG ;
	}
	cld->point_data_size = size;
	cld->point_data = (void *) malloc(cld->nb_points * size);
	if (cld->point_data == NULL) {
		CYAN_ERROR( ERR_MALLOC ) ;
        return ERR_MALLOC ;
	}
	ret = ERR_OK ;
    for (i = 0; i < cld->nb_points; i++) {
		ret += fill_fct( cld, i, context ) ;
	}
	return ret ;
}

int cloud_import_data( cloud_t* cld, size_t data_size, void* data ) {
	if (cld == (cloud_t *) NULL) {
        CYAN_ERROR( ERR_NULL ) ;
		return ERR_NULL ;
	}

	if (cld->point_data == (void *) NULL) {
        cld->point_data = malloc( cld->nb_points * data_size ) ;
        cld->point_data_size = data_size ;

    } else {
        if ( cld->point_data_size != data_size ) {
            free( cld->point_data ) ;
            cld->point_data = malloc( cld->nb_points * data_size ) ;
            cld->point_data_size = data_size ;
        }
    }

    memcpy( cld->point_data, data, cld->nb_points * data_size ) ;
    return ERR_OK ;
}


int cloud_get_data_pointer(cloud_t * cld, int i, void** data_ptr) {
	if ((i<0)||(i>=cld->nb_points)) {
        CYAN_ERROR( ERR_INVALID_ARG ) ;
        return ERR_INVALID_ARG ;
    }
    *data_ptr=cld->point_data+i*cld->point_data_size;
    return ERR_OK ;
}

int cloud_clone( cloud_t* cld, cloud_t** clone ) {
    int i ;

    if ( cld == ( cloud_t* ) NULL ) {
        CYAN_ERROR( ERR_INVALID_ARG ) ;
        return ERR_INVALID_ARG ;
    }

    *clone = cloud_new( cld->nb_points ) ;
    if ( *clone == NULL ) {
        CYAN_ERROR( ERR_MALLOC ) ;
        return ERR_MALLOC ;
    }
    
    memcpy( (*clone)->coords, cld->coords, cld->nb_points*2*sizeof(double) ) ;

    if ( cld->color_X != NULL )
        memcpy( (*clone)->color_X, cld->color_X, cld->nb_points*sizeof(double) ) ;

    if ( cld->color_Y != NULL )
        memcpy( (*clone)->color_Y, cld->color_Y, cld->nb_points*sizeof(double) ) ;
 
    if ( cld->color_Z != NULL )
        memcpy( (*clone)->color_Z, cld->color_Z, cld->nb_points*sizeof(double) ) ;

    (*clone)->illuminant = cld->illuminant ;
    (*clone)->point_data_size = cld->point_data_size ;

    if ( cld->point_data != NULL )
        memcpy( (*clone)->color_Z, cld->color_Z, cld->nb_points*cld->point_data_size ) ;
   

    return ERR_OK ;
}


int cloud_merge( cloud_t* target, cloud_t* cld ) {

    int i;
    int old_size ;
    int new_size ;
    old_size = target->nb_points ;
    new_size = target->nb_points + cld->nb_points ;

    // Recopie des coordonnées

    target->nb_points = new_size ;
    target->coords = realloc( target->coords, new_size*3*sizeof(double) ) ;
    if ( target->coords == NULL ) {
        CYAN_ERROR( ERR_MALLOC ) ;
        return ERR_MALLOC ;
    }
    for ( i=0; i<3*cld->nb_points;i++)
        target->coords[3*old_size+i] = cld->coords[i] ;

    // recopie des couleurs, si possible
    // Sinon, on strip

    if ( (target->color_X != NULL ) &&
            (target->color_Y != NULL ) &&
            (target->color_Z != NULL ) &&
            ( cld->color_X != NULL ) &&
            ( cld->color_Y != NULL ) &&
            ( cld->color_Z != NULL ) &&
            ( target->illuminant == cld->illuminant ) ) {

        target->color_X = realloc( target->color_X, new_size*sizeof(double) ) ;
        if ( target->color_X == NULL ) {
            CYAN_ERROR( ERR_MALLOC ) ;
            return ERR_MALLOC ;
        }
        
        target->color_Y = realloc( target->color_Y, new_size*sizeof(double) ) ;
        if ( target->color_Y == NULL ) {
            CYAN_ERROR( ERR_MALLOC ) ;
            return ERR_MALLOC ;
        }

        target->color_Z = realloc( target->color_Z, new_size*sizeof(double) ) ;
        if ( target->color_Z == NULL ) {
            CYAN_ERROR( ERR_MALLOC ) ;
            return ERR_MALLOC ;
        }

        for ( i=0; i<cld->nb_points;i++) {
            target->color_X[old_size+i] = cld->color_X[i] ;
            target->color_Y[old_size+i] = cld->color_Y[i] ;
            target->color_Z[old_size+i] = cld->color_Z[i] ;
        }

    } else {

        if ( target->color_X != NULL )
            free( target->color_X ) ;
        if ( target->color_Y != NULL )
            free( target->color_Y ) ;
        if ( target->color_Z != NULL )
            free( target->color_Z ) ;
        target->color_X = NULL ;
        target->color_Y = NULL ;
        target->color_Z = NULL ;
    }

    // Recopie des data si possible, sinon on strip
    
   if ( (target->point_data_size == cld->point_data_size) &&
        (target->point_data != NULL) &&
        (cld->point_data != NULL) ) {

    
       target->point_data = realloc( target->point_data, new_size * target->point_data_size ) ;
       if ( target->point_data == NULL ) {
            CYAN_ERROR( ERR_MALLOC ) ;
            return ERR_MALLOC ;
       }    

       memcpy ( target->point_data + old_size*target->point_data_size, cld->point_data, cld->nb_points*cld->point_data_size ) ;

   } else {

        if ( target->point_data != NULL )
            free( target->point_data ) ;
        target->point_data = NULL ;
        target->point_data_size = (size_t) 0 ;

   }

   return ERR_OK ;
}

int cloud_add_point( cloud_t* cloud, double x, double y, double z, double color_x, double color_y, double color_z, void* data ) {

    cloud->coords = realloc( cloud->coords, (cloud->nb_points+1)*3*sizeof(double) ) ;
    if ( cloud->coords == NULL ) {
        CYAN_ERROR( ERR_MALLOC ) ;
        return ERR_MALLOC ;
    }
    cloud->coords[cloud->nb_points] = x ;
    cloud->coords[cloud->nb_points+1] = y ;
    cloud->coords[cloud->nb_points+2] = z ;
    cloud->nb_points += 1 ;
    
    if ( cloud->color_X != NULL ) {
        cloud->color_X = realloc( cloud->color_X, cloud->nb_points*sizeof(double) ) ;
        if ( cloud->color_X == NULL ) {
            CYAN_ERROR( ERR_MALLOC ) ;
            return ERR_MALLOC ;
        }
    }
    cloud->color_X[cloud->nb_points-1] = color_x ;

    if ( cloud->color_Y != NULL ) {
        cloud->color_Y = realloc( cloud->color_Y, cloud->nb_points*sizeof(double) ) ;
        if ( cloud->color_Y == NULL ) {
            CYAN_ERROR( ERR_MALLOC ) ;
            return ERR_MALLOC ;
        }
        cloud->color_Y[cloud->nb_points-1] = color_y ;
    }

    if ( cloud->color_Z != NULL ) {
        cloud->color_Z = realloc( cloud->color_Z, cloud->nb_points*sizeof(double) ) ;
        if ( cloud->color_Z == NULL ) {
            CYAN_ERROR( ERR_MALLOC ) ;
            return ERR_MALLOC ;
        }
        cloud->color_Z[cloud->nb_points-1] = color_z ;
    }

    if ( cloud->point_data != NULL ) {
        if ( data == NULL ) {
            CYAN_ERROR( ERR_INVALID_ARG ) ;
            return ERR_INVALID_ARG ;
        }
        cloud->point_data = realloc( cloud->point_data, cloud->nb_points*cloud->point_data_size ) ;
        if ( cloud->point_data == NULL ) {
            CYAN_ERROR( ERR_MALLOC ) ;
            return ERR_MALLOC ;
        }
        memcpy( cloud->point_data+(cloud->nb_points-1)*cloud->point_data_size, data, cloud->point_data_size ) ;
    }

    return ERR_OK ;
}

void cloud_strip_color( cloud_t* cld ) {
    if ( cld->color_X != NULL ) {
        free( cld->color_X );
    }
    cld->color_X = NULL ;

    if ( cld->color_Y != NULL ) {
        free( cld->color_Y ) ;
    }
    cld->color_Y = NULL ;

    if ( cld->color_Z != NULL ) {
        free( cld->color_Z ) ;
    }
    cld->color_Z = NULL ;
}

void cloud_strip_data( cloud_t* cld ) {
    if ( cld->point_data != NULL )
        free( cld->point_data ) ;
    cld->point_data = NULL ;
    cld->point_data_size = (size_t) 0 ;
}

int cloud_del_point( cloud_t* cloud, int point ) {
    
    if ( ( point < 0 ) || ( point >= cloud->nb_points ) ) {
            CYAN_ERROR( ERR_MALLOC ) ;
            return ERR_MALLOC ;
    }

    cloud->coords[3*point] = cloud->coords[3*(cloud->nb_points-1)] ;
    cloud->coords[3*point+1] = cloud->coords[3*(cloud->nb_points-1)+1] ;
    cloud->coords[3*point+2] = cloud->coords[3*(cloud->nb_points-1)+2] ;
   
    cloud->coords=realloc( cloud->coords, 3*(cloud->nb_points-1)*sizeof(double) ) ;
    if ( cloud->coords == NULL ) {
        CYAN_ERROR( ERR_MALLOC ) ;
        return ERR_MALLOC ;
    }

    if ( cloud->color_X != NULL ) {
        cloud->color_X[point] = cloud->color_X[cloud->nb_points-1] ;
        cloud->color_X=realloc( cloud->color_X, (cloud->nb_points-1)*sizeof(double) ) ;
        if ( cloud->color_X == NULL ) {
            CYAN_ERROR( ERR_MALLOC ) ;
            return ERR_MALLOC ;
        }
    }

    if ( cloud->color_Y != NULL ) {
        cloud->color_Y[point] = cloud->color_Y[cloud->nb_points-1] ;
        cloud->color_Y=realloc( cloud->color_Y, (cloud->nb_points-1)*sizeof(double) ) ;
        if ( cloud->color_Y == NULL ) {
            CYAN_ERROR( ERR_MALLOC ) ;
            return ERR_MALLOC ;
        }
    }
    
    if ( cloud->color_Z != NULL ) {
        cloud->color_Z[point] = cloud->color_Z[cloud->nb_points-1] ;
        cloud->color_Z=realloc( cloud->color_Z, (cloud->nb_points-1)*sizeof(double) ) ;
        if ( cloud->color_Z == NULL ) {
            CYAN_ERROR( ERR_MALLOC ) ;
            return ERR_MALLOC ;
        }
    }

    if( cloud->point_data != NULL) {
        memcpy(cloud->point_data+point*cloud->point_data_size,
               cloud->point_data+(cloud->nb_points-1)*cloud->point_data_size,
               cloud->point_data_size) ;
        cloud->point_data=realloc( cloud->point_data, (cloud->nb_points-1)*cloud->point_data_size ) ;
        if ( cloud->point_data == NULL ) {
            CYAN_ERROR( ERR_MALLOC ) ;
            return ERR_MALLOC ;
        }
    }

   return ERR_OK ;
}





