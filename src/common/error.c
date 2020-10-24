#include "cyan/common/error.h"

#if HAVE_GETTEXT
    
    #include <libintl.h>
    #define _(STRING) dgettext ( "cyan", STRING )

#else

    #define _(STRING) STRING

#endif


// TODO:
//
// Mettre dans la fonction init()
// bindtextdomain ("cyan", "/usr/share/locale";


const char* error_get_text ( int error_code ) {
    switch ( error_code ) {
        case ERR_OK:
            return _("OK") ;
        case ERR_NOPE:
            return _("NOPE") ;
        case ERR_MALLOC:
            return _("Memory allocation failed") ;    
        case ERR_NULL:
            return _("NULL pointer was given") ;
        case ERR_INVALID_ARG:
            return _("Given argument is invalid or has wrong value") ;
        case ERR_NOT_IMPLEMENTED:
            return _("Sorry ... this is not implemented yet :(") ;
        case ERR_FILE:
            return _("Could not open/read or write file. Oops.") ;
        case ERR_FULL:
            return _("Structure is full") ;
        case ERR_EMPTY: 
            return _("Structure is empty") ;
        case ERR_PLUGIN: 
            return _("Plugin could not be loaded") ;
    }
    return _("Unknown error") ;
}

