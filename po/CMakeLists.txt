if(ENABLE_GETTEXT)

    file(STRINGS LINGUAS linguas)

    foreach(lang ${linguas})
        install(FILES ${PROJECT_SOURCE_DIR}/po/${lang}/cyan.mo
            DESTINATION "share/locale/${lang}/LC_MESSAGES/" )
    endforeach()
    
    add_custom_target( create_pot xgettext --keyword=_ 
        --language=C --add-comments  --sort-output 
        -o ${PROJECT_SOURCE_DIR}/po/cyan.pot ${PROJECT_SOURCE_DIR}/src/common/error.c )

    add_custom_target( update_pot xgettext -j --keyword=_ 
        --language=C --add-comments  --sort-output 
        -o ${PROJECT_SOURCE_DIR}/po/cyan.pot ${PROJECT_SOURCE_DIR}/src/common/error.c )

    foreach(lang ${linguas})
        file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/po/${lang})
        add_custom_target( create_po_${lang} 
            msginit --input=${PROJECT_SOURCE_DIR}/po/cyan.pot
            --locale=${lang} --output=${PROJECT_SOURCE_DIR}/po/${lang}/cyan.po )
        add_custom_target( update_po_${lang} 
            msgmerge --update 
            ${PROJECT_SOURCE_DIR}/po/${lang}/cyan.po
            ${PROJECT_SOURCE_DIR}/po/cyan.pot)
        add_custom_target( create_mo_${lang} 
            msgfmt --output-file=${PROJECT_SOURCE_DIR}/po/${lang}/cyan.mo 
            ${PROJECT_SOURCE_DIR}/po/${lang}/cyan.po)
    endforeach()
endif()
