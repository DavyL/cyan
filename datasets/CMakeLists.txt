file( GLOB images "images/*"  )


install(FILES ${images} DESTINATION "${INSTALL_DATASETS_PATH}/images" )  
