#ifndef CYAN_COLOR_H
#define CYAN_COLOR_H

#include <stddef.h>
typedef struct {
	double  X ;                         // XYZ as defined in CIE 1976
	double  Y ;
	double  Z ;

	int monochrome ;                    // if set to zero (remark : isn't it the opposite ?, Leo) , only the Y plane will be used
				
} color_t ;
//---------------------------------------------------------
int color_assign( color_t * dst, double X, double Y, double Z, int monochrome ); 

//---------------------------------------------------------

enum cyan_colorspace {
    CYAN_XYZ        =0,
    CYAN_xyY        =1,
    CYAN_Lab        =2,
    CYAN_LCHab      =3,
    CYAN_Luv        =4,
    CYAN_LCHuv      =5,
    CYAN_RGB        =6,
    CYAN_spectrum   =7
} ;

enum cyan_refwhite {
    CYAN_A      =0,
    CYAN_B      =1,
    CYAN_C      =2,
    CYAN_D50    =3,
    CYAN_D55    =4,
    CYAN_D65    =5,
    CYAN_D75    =6,
    CYAN_E      =7,
    CYAN_F2     =8,
    CYAN_F7     =9,
    CYAN_F11    =10,
} ;

enum cyan_rgb_workspace {
    CYAN_RGB_ADOBE_1998      =0,
    CYAN_RGB_APPLE           =1,
    CYAN_RGB_BEST            =2,
    CYAN_RGB_BETA            =3,
    CYAN_RGB_BRUCE           =4,
    CYAN_RGB_CIE             =5,
    CYAN_RGB_COLORMATCH      =6,
    CYAN_RGB_DON_4           =7,
    CYAN_RGB_ECI             =8,
    CYAN_RGB_EKTA_SPACE_PS5  =9,
    CYAN_RGB_NTSC            =10,
    CYAN_RGB_PAL_SECAM       =11,
    CYAN_RGB_PROPHOTO        =12,
    CYAN_RGB_SMPTE_C         =13,
    CYAN_RGB_sRGB            =14,
    CYAN_RGB_WIDE_GAMUT      =15
} ;

// ---------------------------------------------------------------------------

// All RGB colorspace conversions are made from/to
// D_50 XYZ space

void RGB_to_XYZ( enum cyan_rgb_workspace RGB_space, 
                    double R, double G, double B, 
                    double* X, double* Y, double* Z ) ; 
void XYZ_to_RGB( enum cyan_rgb_workspace RGB_space, 
                    double X, double Y, double Z, 
                    double* R, double* G, double* B ) ;  

void xyY_to_XYZ( double xyY_x, double xyY_y, double xyY_Y, 
                    double* X, double* Y, double* Z ) ;

void XYZ_to_xyY( double X, double Y, double Z,
                    double* xyY_x, double* xyY_y, double* xyY_Y) ; 

void Lab_to_XYZ( enum cyan_refwhite white, 
                    double L, double a, double b, 
                    double* X, double* Y, double* Z) ;
void XYZ_to_Lab( enum cyan_refwhite white, 
                    double X, double Y, double Z,
                    double* L, double* a, double* b ) ;

void LCHab_to_XYZ( enum cyan_refwhite white, 
                    double L, double C, double H, 
                    double* X, double* Y, double* Z) ;
void XYZ_to_LCHab( enum cyan_refwhite white, 
                    double X, double Y, double Z,
                    double* L, double* C, double* H ) ;

void Luv_to_XYZ( enum cyan_refwhite white, 
                    double L, double u, double v, 
                    double* X, double* Y, double* Z) ;
void XYZ_to_Luv( enum cyan_refwhite white, 
                    double X, double Y, double Z,
                    double* L, double* u, double* v ) ;

void LCHuv_to_XYZ( enum cyan_refwhite white, 
                    double L, double C, double H, 
                    double* X, double* Y, double* Z) ;
void XYZ_to_LCHuv( enum cyan_refwhite white, 
                    double X, double Y, double Z,
                    double* L, double* C, double* H ) ;

void spectrum_to_XYZ() ;    // TODO

// ---------------------------------------------------------------------------

void refwhite_get_xyz( enum cyan_refwhite w, double* X, double* Y,double* Z) ;

void rgb_get_name( enum cyan_rgb_workspace rgb, char* str, size_t str_size ) ;
void refwhite_get_name( enum cyan_refwhite w  , char* str, size_t str_size ) ;

#endif
