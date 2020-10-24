#include <math.h>
#include "cyan/color/color.h"

// Watch out for black, where X=Y=Z=0. In that case, you may want to set x and y 
// to the chromaticity coordinates of your reference white.

void XYZ_to_xyY( double X, double Y, double Z,
                    double* xyY_x, double* xyY_y, double* xyY_Y) {

    *xyY_x = X / ( X+Y+Z ) ;
    
    *xyY_y = Y / ( X+Y+Z ) ;
    
    *xyY_Y = Y ;

}
