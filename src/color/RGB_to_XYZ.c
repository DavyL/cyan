#include <math.h>
#include "cyan/color/color.h"

void RGB_to_XYZ(
        enum cyan_rgb_workspace RGB_space,
        double R, double G, double B,
        double* X, double* Y,double* Z) {

    // R,G,B must be in the nominal range [0,1]
    // X,Y,Z will also be in the range[0,1]
    // Ref white for XYZ is D_50 illuminant

    double r,g,b ;  // RGB components made linear to energy
    double k ;      // Used in L* Gamma

    switch( RGB_space ) {
        case CYAN_RGB_ADOBE_1998:
            // Gamma : 2.2
            // Ref White: D65
            // Adaptation Method : Bradford
            r = pow( R, 2.2 ) ;
            g = pow( G, 2.2 ) ;
            b = pow( B, 2.2 ) ;
            *X = 0.6097559*r+0.2052401*g+0.1492240*b;
            *Y = 0.3111242*r+0.6256560*g+0.0632197*b;
            *Z = 0.0194811*r+0.0608902*g+0.7448387*b;
            break ;
        case CYAN_RGB_APPLE:
            // Gamma : 1.8
            // Ref White: D65
            // Adaptation Method : Bradford
            r = pow( R, 1.8 ) ;
            g = pow( G, 1.8 ) ;
            b = pow( B, 1.8 ) ;
            *X=0.4755678*r+0.3396722*g+0.1489800*b;
            *Y=0.2551812*r+0.6725693*g+0.0722496*b;
            *Z=0.0184697*r+0.1133771*g+0.6933632*b;
            break ;
        case CYAN_RGB_BEST:
            // Gamma : 2.2
            // Ref White: D50
            // Adaptation Method : None
            r = pow( R, 2.2 ) ;
            g = pow( G, 2.2 ) ;
            b = pow( B, 2.2 ) ;
            *X=0.6326696*r+0.2045558*g+0.1269946*b;
            *Y=0.2284569*r+0.7373523*g+0.0341908*b;
            *Z=0.0000000*r+0.0095142*g+0.8156958*b;
            break ;
        case CYAN_RGB_BETA:
            // Gamma : 2.2
            // Ref White: D50
            // Adaptation Method : None
            r = pow( R, 2.2 ) ;
            g = pow( G, 2.2 ) ;
            b = pow( B, 2.2 ) ;
            *X=0.6712537*r+0.1745834*g+0.1183829*b;
            *Y=0.3032726*r+0.6637861*g+0.0329413*b;
            *Z=0.0000000*r+0.0407010*g+0.7845090*b;
            break ;
        case CYAN_RGB_BRUCE:
            // Gamma : 2.2
            // Ref White: D65 
            // Adaptation Method : Bradford
            r = pow( R, 2.2 ) ;
            g = pow( G, 2.2 ) ;
            b = pow( B, 2.2 ) ;
            *X=0.4941816*r+0.3204834*g+0.1495550*b;
            *Y=0.2521531*r+0.6844869*g+0.0633600*b;
            *Z=0.0157886*r+0.0629304*g+0.7464909*b;
            break ;
        case CYAN_RGB_CIE:
            // Gamma : 2.2
            // Ref White: E
            // Adaptation Method : Bradford
            r = pow( R, 2.2 ) ;
            g = pow( G, 2.2 ) ;
            b = pow( B, 2.2 ) ;
            *X=0.4868870*r+0.3062984*g+0.1710347*b;
            *Y=0.1746583*r+0.8247541*g+0.0005877*b;
            *Z=-0.0012563*r+0.0169832*g+0.8094831*b;
            break ;
        case CYAN_RGB_COLORMATCH:
            // Gamma : 1.8
            // Ref White: D50
            // Adaptation Method : None
            r = pow( R, 1.8 ) ;
            g = pow( G, 1.8 ) ;
            b = pow( B, 1.8 ) ;
            *X=0.5093439*r+0.3209071*g+0.1339691*b;
            *Y=0.2748840*r+0.6581315*g+0.0669845*b;
            *Z=0.0242545*r+0.1087821*g+0.6921735*b;
            break ;
        case CYAN_RGB_DON_4:
            // Gamma : 2.2
            // Ref White: D50
            // Adaptation Method : None
            r = pow( R, 2.2 ) ;
            g = pow( G, 2.2 ) ;
            b = pow( B, 2.2 ) ;
            *X=0.6457711*r+0.1933511*g+0.1250978*b;
            *Y=0.2783496*r+0.6879702*g+0.0336802*b;
            *Z=0.0037113*r+0.0179861*g+0.8035125*b;
            break ;
        case CYAN_RGB_ECI:
            // Gamma : L*
            // Ref White: D50
            // Adaptation Method : None
            k = 903.3 ; // As defined in CIE standard
            if (R<= 0.08) {
                r = 100.0 * R / k ;
            } else {
                r = pow(((R+0.16)/1.16),3.0) ;
            }
            if (G<= 0.08) {
                g = 100.0 * G / k ;
            } else {
                g = pow(((G+0.16)/1.16),3.0) ;
            }
            if (B<= 0.08) {
                b = 100.0 * B / k ;
            } else {
                b = pow(((B+0.16)/1.16),3.0) ;
            }
            *X=0.6502043*r+0.1780774*g+0.1359384*b;
            *Y=0.3202499*r+0.6020711*g+0.0776791*b;
            *Z=0.0000000*r+0.0678390*g+0.7573710*b;
            break;
        case CYAN_RGB_EKTA_SPACE_PS5:
            // Gamma : 2.2
            // Ref White: D50
            // Adaptation Method : None
            r = pow( R, 2.2 ) ;
            g = pow( G, 2.2 ) ;
            b = pow( B, 2.2 ) ;
            *X=0.5938914*r+0.2729801*g+0.0973485*b;
            *Y=0.2606286*r+0.7349465*g+0.0044249*b;
            *Z=0.0000000*r+0.0419969*g+0.7832131*b;
            break;
        case CYAN_RGB_NTSC:
            // Gamma : 2.2
            // Ref White: C
            // Adaptation Method : Bradford
            r = pow( R, 2.2 ) ;
            g = pow( G, 2.2 ) ;
            b = pow( B, 2.2 ) ;
            *X=0.6343706*r+0.1852204*g+0.1446290*b;
            *Y=0.3109496*r+0.5915984*g+0.0974520*b;
            *Z=-0.0011817*r+0.0555518*g+0.7708399*b;
            break;
        case CYAN_RGB_PAL_SECAM:
            // Gamma : 2.2
            // Ref White: D65
            // Adaptation Method : Bradford
            r = pow( R, 2.2 ) ;
            g = pow( G, 2.2 ) ;
            b = pow( B, 2.2 ) ;
            *X=0.4552773*r+0.3675500*g+0.1413926*b;
            *Y=0.2323025*r+0.7077956*g+0.0599019*b;
            *Z=0.0145457*r+0.1049154*g+0.7057489*b;
            break;
        case CYAN_RGB_PROPHOTO:
            // Gamma : 1.8
            // Ref White: D50
            // Adaptation Method : None
            r = pow( R, 1.8 ) ;
            g = pow( G, 1.8 ) ;
            b = pow( B, 1.8 ) ;
            *X=0.7976749*r+0.1351917*g+0.0313534*b;
            *Y=0.2880402*r+0.7118741*g+0.0000857*b;
            *Z=0.0000000*r+0.0000000*g+0.8252100*b;
            break;
        case CYAN_RGB_SMPTE_C:
            // Gamma : 2.2
            // Ref White: D65
            // Adaptation Method : Bradford
            r = pow( R, 2.2 ) ;
            g = pow( G, 2.2 ) ;
            b = pow( B, 2.2 ) ;
            *X=0.4163290*r+0.3931464*g+0.1547446*b;
            *Y=0.2216999*r+0.7032549*g+0.0750452*b;
            *Z=0.0136576*r+0.0913604*g+0.7201920*b;
            break;
        case CYAN_RGB_sRGB:
            // Gamma : sRGB (approx. 2.2)
            // Ref White: D65
            // Adaptation Method : Bradford
            if (R<= 0.04045) {
                r = R / 12.92 ;
            } else {
                r = pow(((R+0.055)/1.055),2.4) ;
            }
            if (G<= 0.04045) {
                g = G / 12.92 ;
            } else {
                g = pow(((G+0.055)/1.055),2.4) ;
            }
            if (B<= 0.04045) {
                b = B / 12.92 ;
            } else {
                b = pow(((B+0.055)/1.055),2.4) ;
            }
            *X=0.4360747*r+0.3850649*g+0.1430804*b;
            *Y=0.2225045*r+0.7168786*g+0.0606169*b;
            *Z=0.0139322*r+0.0971045*g+0.7141733*b;
            break;
        case CYAN_RGB_WIDE_GAMUT:
            // Gamma : 2.2
            // Ref White: D50
            // Adaptation Method : None
            r = pow( R, 2.2 ) ;
            g = pow( G, 2.2 ) ;
            b = pow( B, 2.2 ) ;
            *X=0.7161046*r+0.1009296*g+0.1471858*b;
            *Y=0.2581874*r+0.7249378*g+0.0168748*b;
            *Z=0.0000000*r+0.0517813*g+0.7734287*b;
            break;
        default:
            break ;
    }
}
