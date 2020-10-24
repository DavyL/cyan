#include <math.h>
#include "cyan/color/color.h"

#define EPSILON 0.008856
#define KAPPA   903.3


void Lab_to_XYZ( enum cyan_refwhite white, 
                    double L, double a, double b, 
                    double* X, double* Y, double* Z) {

    double Xr, Yr, Zr ;
    double fx, fy, fz ;
    double xr, yr, zr ; 
    refwhite_get_xyz( white, &Xr, &Yr, &Zr) ;

    fy = ( L + 16.0 ) / 116.0 ;
    fx = a / 500.0 + fy ;
    fz = fy - b / 200.0 ;

   if ( pow(fx,3.0) > EPSILON )
       xr = pow(fx,3.0) ;
   else
       xr = (116*fx-16)/KAPPA ;

   if ( L > KAPPA*EPSILON )
       yr = pow( (L+16.0)/116.0, 3.0 ) ;
   else
       yr = L / KAPPA ;

   if ( pow(fz,3.0) > EPSILON )
       zr = pow(fz, 3.0 ) ;
   else
       zr = ( 116.0 * fz - 16.0 ) / KAPPA ;

    *X = xr * Xr ;
    *Y = yr * Yr ;
    *Z = zr * Zr ;

}
