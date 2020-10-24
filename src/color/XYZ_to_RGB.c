#include <math.h>
#include "cyan/color/color.h"

void XYZ_to_RGB( enum cyan_rgb_workspace RGB_space, 
                    double X, double Y, double Z, 
                    double* R, double* G, double* B ) {

    double r ;
    double g ;
    double b ;

    switch( RGB_space ) {
        case CYAN_RGB_ADOBE_1998:
            // Gamma : 2.2
            // Ref White: D65
            // Adaptation Method : Bradford
            r = 1.9624274*X-0.6105343*Y-0.3413404*Z ; 
            g =-0.9787684*X+1.9161415*Y+0.0334540*Z ;
            b = 0.0286869*X-0.1406752*Y+1.3487655*Z ;
            *R = pow( r, 1.0 / 2.2 ) ;
            *G = pow( g, 1.0 / 2.2 ) ;
            *B = pow( b, 1.0 / 2.2 ) ;
            break ;
        case CYAN_RGB_APPLE:
            // Gamma : 1.8
            // Ref White: D65
            // Adaptation Method : Bradford
            r = 2.8510695*X-1.3605261*Y-0.4708281*Z;
            g =-1.0927680*X+2.0348871*Y+0.0227598*Z;
            b = 0.1027403*X-0.2964984*Y+1.4510659*Z;
            *R = pow( r, 1.0 / 1.8 ) ;
            *G = pow( g, 1.0 / 1.8 ) ;
            *B = pow( b, 1.0 / 1.8 ) ;
            break ;
        case CYAN_RGB_BEST:
            // Gamma : 2.2
            // Ref White: D50
            // Adaptation Method : None
            r = 1.7552599*X-0.4836786*Y-0.2530000*Z ; 
            g =-0.5441336*X+1.5068789*Y+0.0215528*Z ;
            b = 0.0063467*X-0.0175761*Y+1.2256959*Z ;
            *R = pow( r, 1.0 / 2.2 ) ;
            *G = pow( g, 1.0 / 2.2 ) ;
            *B = pow( b, 1.0 / 2.2 ) ;
            break ;
        case CYAN_RGB_BETA:
            // Gamma : 2.2
            // Ref White: D50
            // Adaptation Method : None
            r =  1.6832270*X-0.4282363*Y-0.2360185*Z ; 
            g = -0.7710229*X+1.7065571*Y+0.0446900*Z ;
            b =  0.0400013*X-0.0885376*Y+1.2723640*Z ;
            *R = pow( r, 1.0 / 2.2 ) ;
            *G = pow( g, 1.0 / 2.2 ) ;
            *B = pow( b, 1.0 / 2.2 ) ;
            break ;
        case CYAN_RGB_BRUCE:
            // Gamma : 2.2
            // Ref White: D65 
            // Adaptation Method : Bradford
            r = 2.6502856*X-1.2014485*Y-0.4289936*Z;
            g =-0.9787684*X+1.9161415*Y+0.0334540*Z;
            b = 0.0264570*X-0.1361227*Y+1.3458542*Z;
            *R = pow( r, 1.0 / 2.2 ) ;
            *G = pow( g, 1.0 / 2.2 ) ;
            *B = pow( b, 1.0 / 2.2 ) ;
            break ;
        case CYAN_RGB_CIE:
            // Gamma : 2.2
            // Ref White: E
            // Adaptation Method : Bradford
            r = 2.3638081*X-0.8676030*Y-0.4988161*Z;
            g =-0.5005940*X+1.3962369*Y+0.1047562*Z;
            b = 0.0141712*X-0.0306400*Y+1.2323842*Z;
            *R = pow( r, 1.0 / 2.2 ) ;
            *G = pow( g, 1.0 / 2.2 ) ;
            *B = pow( b, 1.0 / 2.2 ) ;
            break ;
        case CYAN_RGB_COLORMATCH:
            // Gamma : 1.8
            // Ref White: D50
            // Adaptation Method : None
            r =  2.6422874*X-1.2234270*Y-0.3930143*Z;
            g = -1.1119763*X+2.0590183*Y+0.0159614*Z;
            b =  0.0821699*X-0.2807254*Y+1.4559877*Z;
            *R = pow( r, 1.0 / 1.8 ) ;
            *G = pow( g, 1.0 / 1.8 ) ;
            *B = pow( b, 1.0 / 1.8 ) ;
            break ;
        case CYAN_RGB_DON_4:
            // Gamma : 2.2
            // Ref White: D50
            // Adaptation Method : None
            r = 1.7603902*X-0.4881198*Y-0.2536126*Z; 
            g =-0.7126288*X+1.6527432*Y+0.0416715*Z;
            b = 0.0078207*X-0.0347411*Y+1.2447743*Z;
            *R = pow( r, 1.0 / 2.2 ) ;
            *G = pow( g, 1.0 / 2.2 ) ;
            *B = pow( b, 1.0 / 2.2 ) ;
            break ;
        case CYAN_RGB_ECI:
            // Gamma : L*
            // Ref White: D50
            // Adaptation Method : None
            r = 1.7827618*X-0.4969847*Y-0.2690101*Z;
            g =-0.9593623*X+1.9477962*Y-0.0275807*Z;
            b = 0.0859317*X-0.1744674*Y+1.3228273*Z;
            if ( r <= 0.008856 )
                *R = r * 903.3 / 100.0 ;
            else
                *R = 1.16*pow(r,1.0/3.0) - 0.16 ;
            if ( g <= 0.008856 )
                *G = g * 903.3 / 100.0 ;
            else
                *G = 1.16*pow(g,1.0/3.0) - 0.16 ;
            if ( b <= 0.008856 )
                *B = b * 903.3 / 100.0 ;
            else
                *B = 1.16*pow(b,1.0/3.0) - 0.16 ;
            break;
        case CYAN_RGB_EKTA_SPACE_PS5:
            // Gamma : 2.2
            // Ref White: D50
            // Adaptation Method : None
            r = 2.0043819*X-0.7304844*Y-0.2450052*Z;
            g =-0.7110285*X+1.6202126*Y+0.0792227*Z;
            b = 0.0381263*X-0.0868780*Y+1.2725438*Z;
            *R = pow( r, 1.0 / 2.2 ) ;
            *G = pow( g, 1.0 / 2.2 ) ;
            *B = pow( b, 1.0 / 2.2 ) ;
            break;
        case CYAN_RGB_NTSC:
            // Gamma : 2.2
            // Ref White: C
            // Adaptation Method : Bradford
            r = 1.8464881*X-0.5521299*Y-0.2766458*Z;
            g =-0.9826630*X+2.0044755*Y-0.0690396*Z;
            b = 0.0736477*X-0.1453020*Y+1.3018376*Z;
            *R = pow( r, 1.0 / 2.2 ) ;
            *G = pow( g, 1.0 / 2.2 ) ;
            *B = pow( b, 1.0 / 2.2 ) ;
            break;
        case CYAN_RGB_PAL_SECAM:
            // Gamma : 2.2
            // Ref White: D65
            // Adaptation Method : Bradford
            r = 2.9603944*X-1.4678519*Y-0.4685105*Z;
            g =-0.9787684*X+1.9161415*Y+0.0334540*Z;
            b = 0.0844874*X-0.2545973*Y+1.4216174*Z;
            *R = pow( r, 1.0 / 2.2 ) ;
            *G = pow( g, 1.0 / 2.2 ) ;
            *B = pow( b, 1.0 / 2.2 ) ;
            break;
        case CYAN_RGB_PROPHOTO:
            // Gamma : 1.8
            // Ref White: D50
            // Adaptation Method : None
            r = 1.3459433*X-0.2556075*Y-0.0511118*Z;
            g =-0.5445989*X+1.5081673*Y+0.0205351*Z;
            b = 0.0000000*X+0.0000000*Y+1.2118128*Z;
            *R = pow( r, 1.0 / 1.8 ) ;
            *G = pow( g, 1.0 / 1.8 ) ;
            *B = pow( b, 1.0 / 1.8 ) ;
            break;
        case CYAN_RGB_SMPTE_C:
            // Gamma : 2.2
            // Ref White: D65
            // Adaptation Method : Bradford
            r = 3.3921940*X-1.8264027*Y-0.5385522*Z;
            g =-1.0770996*X+2.0213975*Y+0.0207989*Z;
            b = 0.0723073*X-0.2217902*Y+1.3960932*Z;
            *R = pow( r, 1.0 / 2.2 ) ;
            *G = pow( g, 1.0 / 2.2 ) ;
            *B = pow( b, 1.0 / 2.2 ) ;
            break;
        case CYAN_RGB_sRGB:
            // Gamma : sRGB (approx. 2.2)
            // Ref White: D65
            // Adaptation Method : Bradford
            r = 3.1338561*X-1.6168667*Y-0.4906146*Z;
            g =-0.9787684*X+1.9161415*Y+0.0334540*Z;
            b = 0.0719453*X-0.2289914*Y+1.4052427*Z;
            if (r<= 0.0031308) {
                *R=12.92*r ;
            } else {
                *R = 1.055*pow(r,1.0/2.4) - 0.055 ;
            }
            if (g<= 0.0031308) {
                *G=12.92*g ;
            } else {
                *G = 1.055*pow(g,1.0/2.4) - 0.055 ;
            }
            if (b<= 0.0031308) {
                *B=12.92*b ;
            } else {
                *B = 1.055*pow(b,1.0/2.4) - 0.055 ;
            } 
            break;
        case CYAN_RGB_WIDE_GAMUT:
            // Gamma : 2.2
            // Ref White: D50
            // Adaptation Method : None
            r = 1.4628067*X-0.1840623*Y-0.2743606*Z;
            g =-0.5217933*X+1.4472381*Y+0.0677227*Z;
            b = 0.0349342*X-0.0968930*Y+1.2884099*Z;
            *R = pow( r, 1.0 / 2.2 ) ;
            *G = pow( g, 1.0 / 2.2 ) ;
            *B = pow( b, 1.0 / 2.2 ) ;
            break;
        default:
            break ;
    }
}


