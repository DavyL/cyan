#include <math.h>
#include <cyan/color/color.h>

#define EPSILON 0.008856
#define KAPPA   903.3

void XYZ_to_LCHab( enum cyan_refwhite white, 
                    double X, double Y, double Z,
                    double* L, double* C, double* H ) {

    double Xr, Yr, Zr ;
    double xr, yr, zr ;
    double fx, fy, fz ;
    double a, b ;

    refwhite_get_xyz( white, &Xr, &Yr, &Zr) ;

    xr = X / Xr ;
    yr = Y / Yr ;
    zr = Z / Zr ;

    if ( xr > EPSILON ) {
        fx = pow( xr, 1.0 / 3.0 ) ;
    } else {
        fx = ( KAPPA * xr + 16.0 ) / 116.0 ;
    }
    
    if ( yr > EPSILON ) {
        fy = pow( yr, 1.0 / 3.0 ) ;
    } else {
        fy = ( KAPPA * yr + 16.0 ) / 116.0 ;
    }

    if ( zr > EPSILON ) {
        fz = pow( zr, 1.0 / 3.0 ) ;
    } else {
        fz = ( KAPPA * zr + 16.0 ) / 116.0 ;
    }

    *L = 116.0 * fy - 16.0 ; 
    a = 500.0 * ( fx - fy ) ;
    b = 200.0 * ( fy - fz ) ;

    *C = sqrt( pow(a, 2.0 ) + pow(b,2.0) ) ;
    
    *H = atan2( a, b ) ;
    *H = *H / M_PI * 180.0 ;
    if ( *H < 0 )
        *H = *H + 360 ;

}
