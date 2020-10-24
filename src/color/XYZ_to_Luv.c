#include <math.h>
#include <cyan/color/color.h>

#define EPSILON 0.008856
#define KAPPA   903.3

void XYZ_to_Luv( enum cyan_refwhite white, 
                    double X, double Y, double Z,
                    double* L, double* u, double* v ) {

    double Xr, Yr, Zr ;
    double ur, vr, yr ;
    double uu, vv ;
    refwhite_get_xyz( white, &Xr, &Yr, &Zr) ;

    ur = 4.0*Xr / ( Xr + 15.0*Yr + 3.0*Zr ) ;
    vr = 9.0*Yr / ( Xr + 15.0*Yr + 3.0*Zr ) ;

    uu = 4.0*X / ( X + 15.0*Y + 3.0*Z ) ;
    vv = 9.0*Y / ( X + 15.0*Y + 3.0*Z ) ;

    yr = Y / Yr ;

    if ( yr > EPSILON ) 
        *L = 116.0 * pow( yr, 1.0/ 3.0 ) - 16.0 ;
    else
        *L = KAPPA * yr ;

    *u = 13.0 * (*L) * ( uu - ur ) ; 
    *v = 13.0 * (*L) * ( vv - vr ) ; 

}
