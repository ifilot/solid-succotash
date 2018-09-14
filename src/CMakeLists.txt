#/**************************************************************************
 #  CMakeLists.txt  --  This file is part of DFTCXX.                     *
 #                                                                       *
 #  Copyright (C) 2016, Ivo Filot                                        *
 #                                                                       *
 #  DFTCXX is free software:                                             *
 #  you can redistribute it and/or modify it under the terms of the      *
 #  GNU General Public License as published by the Free Software         *
 #  Foundation, either version 3 of the License, or (at your option)     *
 #  any later version.                                                   *
 #                                                                       *
 #  DFTCXX is distributed in the hope that it will be useful,            *
 #  but WITHOUT ANY WARRANTY; without even the implied warranty          *
 #  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 #  See the GNU General Public License for more details.                 *
 #                                                                       *
 #  You should have received a copy of the GNU General Public License    *
 #  along with this program.  If not, see http://www.gnu.org/licenses/.  *
 #                                                                       *
#**************************************************************************/

# set minimum cmake requirements
cmake_minimum_required(VERSION 2.8)
project (sme)

# add custom directory to look for .cmake files
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules )

# Enable release build
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

find_package(PkgConfig REQUIRED)
pkg_check_modules(CAIRO cairo REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_BINARY_DIR}
                    ${CAIRO_INCLUDE_DIRS})

# Add sources
file(GLOB SOURCES "*.cpp")

# Set executable
add_executable(sme ${SOURCES})

# link libraries
target_link_libraries(sme cairo)

###
# Installing
##
install (TARGETS sme DESTINATION bin)
