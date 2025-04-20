cmake_minimum_required(VERSION 3.16)
project(T001_geo)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

# モジュールを有効化する
igl_include(core)

# ソースファイルを集める 
file(GLOB SRC_FILES *.cpp)

# 実行ファイルを生成する 
add_executable(${PROJECT_NAME} ${SRC_FILES})

#
add_subdirectory(lib gmsh)

# 必要なライブラリをリンクする 
target_link_libraries(${PROJECT_NAME} PUBLIC 
  igl::core
  gmsh
)

# データへのパスを指定する 
target_compile_definitions( ${PROJECT_NAME} PRIVATE DATA_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data\" )

