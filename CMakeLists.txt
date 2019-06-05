cmake_minimum_required(VERSION 3.5 FATAL_ERROR)  # 3.5 avoid invalid C compiler

project(mp2d
        VERSION 1.0
        LANGUAGES CXX)
set(mp2d_AUTHORS      "Chandler S. Greenwell and Gregory J. O. Beran")
set(mp2d_DESCRIPTION  "Dispersion correction to MP2")
set(mp2d_URL          "https://github.com/Chandemonium/MP2D.git")
set(mp2d_LICENSE      "MIT")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

################  Options: Overview and Not Otherwise Mentioned  ###############

#  <<<  CMake build overview  >>>
#
#    >>> ls
#    MP2D.h LICENSE mp2d_paramseters.dat ...
#    >>> cmake -H. -Bobjdir -DCMAKE_INSTALL_PREFIX=/path/to/install-mp2d ...
#    ...
#    -- Generating done
#    -- Build files have been written to: /current/dir/objdir
#    >>> cd objdir && make -j`getconf _NPROCESSORS_ONLN`
#    >>> make install

#  <<<  Required build dependencies that MP2D can't build itself  >>>
#
#    - CMake (e.g., `conda install cmake`)
#    - C++ compiler

############################  Options: Build How?  #############################
include(psi4OptionsTools)

option_with_default(CMAKE_BUILD_TYPE "Build type (Release or Debug)" Release)

#################################  Main Project  ################################
include(GNUInstallDirs)
include(autocmake_safeguards)

# <<<  Build >>>

configure_file(MP2D.cpp
               ${CMAKE_CURRENT_BINARY_DIR}/MP2D.cpp @ONLY)
configure_file(mp2d_parameters.dat
               ${CMAKE_CURRENT_BINARY_DIR}/mp2d_parameters.dat @ONLY)

set(sources_list Main.cpp
                 ${CMAKE_CURRENT_BINARY_DIR}/MP2D.cpp)

add_executable(mp2d ${sources_list})
target_include_directories(mp2d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# <<<  Install  >>>

install(TARGETS mp2d
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES CKS_C6coeffs.dat 
              UCHF_C6coeffs.dat
              ${CMAKE_CURRENT_BINARY_DIR}/mp2d_parameters.dat 
        DESTINATION ${CMAKE_INSTALL_DATADIR}/mp2d)
