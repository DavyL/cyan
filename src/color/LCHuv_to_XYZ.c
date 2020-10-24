#include <math.h>
#include "cyan/color/color.h"

#define EPSILON 0.008856
#define KAPPA   903.3

// h est exprimé en degrés

void LCHuv_to_XYZ( enum cyan_refwhite white, 
                    double L, double C, double H, 
                    double* X, double* Y, double* Z) {

    double Xr, Yr, Zr ;
    double a,b,c,d ;
    double u0, v0 ;
    double u ;
    double v ;

    u = C*cos(H * M_PI / 180.0 ) ;
    v = C*sin(H * M_PI / 180.0 ) ;
   

    refwhite_get_xyz( white, &Xr, &Yr, &Zr) ;

    u0 = 4.0*Xr/( Xr + 15.0*Yr + 3.0*Zr ) ;
    v0 = 9.0*Yr/( Xr + 15.0*Yr + 3.0*Zr ) ;


   if ( L > KAPPA*EPSILON )
       *Y = pow( (L+16.0)/116.0, 3.0 ) ;
   else
       *Y = L / KAPPA ;

    a = 1.0/3.0*( 52.0 * L / ( u + 13.0*L*u0 ) - 1.0 ) ;
    b = -5.0*(*Y) ;
    c = -1.0 / 3.0 ;
    d = (*Y) * ( 39.0* L / ( v + 13.0 * L * v0 ) - 5.0 ) ;
   
    *X = ( d - b ) / ( a - c ) ;
    *Z = (*X)*a + b ;

}
