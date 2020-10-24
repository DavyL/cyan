#include <math.h>
#include "cyan/color/color.h"

// The output (X,Y,Z) values are in the same nominal range as the input xyY_Y


void xyY_to_XYZ( double xyY_x, double xyY_y, double xyY_Y, double* X, double* Y, double* Z ) {

    if ( xyY_y == 0.0 ) {
        *X = 0.0 ;
        *Y = 0.0 ;
        *Z = 0.0 ;
    } else {
        *X = xyY_x * xyY_Y / xyY_y ;
        *Y = xyY_Y ;
        *Z = ( 1 - xyY_x - xyY_y ) * xyY_Y / xyY_y ;
    }
}
