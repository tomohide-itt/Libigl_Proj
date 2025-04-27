# トップレベルプロジェクトかどうかを検出
get_directory_property(LIBIGL_PARENT_DIR PARENT_DIRECTORY)
if(NOT LIBIGL_PARENT_DIR)
    set(LIBIGL_TOPLEVEL_PROJECT ON)
else()
    set(LIBIGL_TOPLEVEL_PROJECT OFF)
endif()

# 必要なCMakeのバージョンをチェック
# トップレベルプロジェクトであれば，CMakeの最低必要なバージョンを設定
# そうでない場合，CMakeのバージョンがREQUIRED_CMAKE_VERSION以下であれば，エラーメッセージを表示して終了
set(REQUIRED_CMAKE_VERSION "3.16.0")
if(LIBIGL_TOPLEVEL_PROJECT)
    cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
else()
    if(${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
        message(FATAL_ERROR "CMake required version to build Libigl is ${REQUIRED_CMAKE_VERSION}")
    endif()
endif()

################################################################################
project(libigl VERSION 2.5.0)

# CMake module path
list(PREPEND CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_LIST_DIR}/cmake/
    ${CMAKE_CURRENT_LIST_DIR}/cmake/igl
    ${CMAKE_CURRENT_LIST_DIR}/cmake/find
    ${CMAKE_CURRENT_LIST_DIR}/cmake/recipes/external
)

set_property(GLOBAL PROPERTY __igl_module_path ${CMAKE_MODULE_PATH})

# Set default options
set(LIBIGL_DEFAULT_CGAL ${LIBIGL_TOPLEVEL_PROJECT})

# Build tutorials
option(LIBIGL_BUILD_TUTORIALS  "Build libigl tutorial"                 ${LIBIGL_TOPLEVEL_PROJECT})

# Permissive modules. These modules are available under MPL2 license, and their dependencies are available
# under a permissive or public domain license.
option(LIBIGL_GLFW             "Build target igl::glfw"             ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_OPENGL           "Build target igl::opengl"           ${LIBIGL_TOPLEVEL_PROJECT})

# Copyleft modules. These modules are available under GPL license, and their dependencies are
# available under a copyleft license.
option(LIBIGL_COPYLEFT_CORE    "Build target igl_copyleft::core"    ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_COPYLEFT_CGAL    "Build target igl_copyleft::cgal"    ${LIBIGL_DEFAULT_CGAL})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")


# Include CMake helper functions
include(igl_add_library)
include(igl_add_test)
include(igl_add_tutorial)
include(igl_copy_dll)
include(igl_include)
include(igl_install)
include(igl_target_sources)

# Defines CMake targets for selected libigl modules
include(libigl)

if(LIBIGL_BUILD_TUTORIALS)
    add_subdirectory(tutorial)
endif()

