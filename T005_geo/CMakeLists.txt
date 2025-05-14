# 必要なCMakeのバージョンをチェック
cmake_minimum_required(VERSION 3.16)

################################################################################
# プロジェクトの名前とバージョンを設定
project(libigl VERSION 2.5.0)
set(PRGNAME "T005_geo")

# CMakeのモジュール検索パス（CMAKE_MODULE_PATH）に複数のディレクトリを先頭に追加（PREPEND）
# これにより，CMakeは指定されたディレクトリ内のモジュールを優先的に検索する
# CMAKE_CURRENT_LIST_DIRは現在のCMakeLists.txtファイルのディレクトリの絶対パスを示す
list(PREPEND CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_LIST_DIR}/cmake/igl
    ${CMAKE_CURRENT_LIST_DIR}/cmake/recipes/external
)

# CMakeのモジュールパスをグローバルプロパティ__igl_module_pathに設定
# これにより，libiglのモジュールを他のCMakeLists.txtファイルからも参照できるようになる
# __igl_module_pathはlibiglのモジュールを検索するためのパスを格納するプロパティ
set_property(GLOBAL PROPERTY __igl_module_path ${CMAKE_MODULE_PATH})

# ディレクトリパスの設定 
# CMAKE_BINARY_DIRは，CMakeがビルドを行うために使用するディレクトリのパスを示す
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/libs" )

# CMake ヘルパー関数をインクルード
# これらのヘルパー関数は，libiglのビルドプロセスを簡素化するために使用される
include(igl_add_library)
include(igl_add_executable)
include(igl_copy_dll)
include(igl_include)
include(igl_install)
include(igl_target_sources)

# libiglのモジュールをCMakeのターゲットとして定義する
igl_include(core)
igl_include(glfw)
igl_include(copyleft cgal)

# 自作ライブラリを追加
add_subdirectory( libs/gmsh )
add_subdirectory( libs/migl )

# 実行可能ファイルの追加
igl_add_executable(${PRGNAME} igl::glfw igl_copyleft::cgal gmsh migl )

target_compile_features(${PRGNAME} PUBLIC cxx_std_17)

target_include_directories( ${PRGNAME} PRIVATE
    ${LIB_DIR}/gmsh
    ${LIB_DIR}/migl
)

# データへのパスを指定する 
target_compile_definitions( ${PRGNAME} PRIVATE DATA_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data\" )