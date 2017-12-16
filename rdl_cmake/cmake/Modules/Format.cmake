## Some fancy cmake macros courtesy of the folks at Houston Mechatronics ##

####
# Macro(s) that build cmake targets that automatigically reformat code when run
#
####

# Uncrustify is a style beautifier for c++ code
find_program(UNCRUSTIFY_EXE uncrustify)

set(run_cpp_format ON)
if(NOT UNCRUSTIFY_EXE)
    message(WARNING "Uncrustify config file not found! Can't autoformat!")
    set(run_cpp_format OFF)
else()
    if(EXISTS "${PROJECT_SOURCE_DIR}/config/uncrustify.cfg")
        set(uncrustify_cfg "${PROJECT_SOURCE_DIR}/config/uncrustify.cfg")
        set(run_cpp_format ON)
        message(STATUS "Project ${PROJECT_NAME} using uncrustify config file in ${PROJECT_SOURCE_DIR}/config/uncrustify/uncrustify.cfg")
    else()
        set(run_cpp_format OFF)
        message(WARNING "No uncrustify config found for ${PROJECT_NAME}!")
    endif()
endif()

macro(format_code_target)

    file(GLOB_RECURSE file_list ${PROJECT_NAME}/src/*.c* include/*/*.h*)

    # If there are no c++ files, don't run uncrustify
    if(NOT file_list)
        set(run_cpp_format OFF)
    endif()

    if(${run_cpp_format})
        add_custom_target(format_${PROJECT_NAME}
                COMMAND "${UNCRUSTIFY_EXE}" --replace --no-backup -l "CPP" -c "${uncrustify_cfg}" ${file_list}
                WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
                COMMENT "Uncrustifying c++ files for ${PROJECT_NAME}")
    else()
        message(WARNING "Unable to format C++ source code for ${PROJECT_NAME}. Check that uncrustify is installed and that you have an uncrustify.cfg file!")
    endif()
endmacro()
