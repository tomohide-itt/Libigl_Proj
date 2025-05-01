# トップレベルプロジェクトかどうかを検出
get_directory_property(LIBIGL_PARENT_DIR PARENT_DIRECTORY)
if(NOT LIBIGL_PARENT_DIR)
    set(LIBIGL_TOPLEVEL_PROJECT ON)
else()
    set(LIBIGL_TOPLEVEL_PROJECT OFF)
endif()

# 必要なCMakeのバージョンをチェック
cmake_minimum_required(VERSION 3.16)

################################################################################
# プロジェクトの名前とバージョンを設定
project(libigl VERSION 2.5.0)

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

# Set default options
set(LIBIGL_DEFAULT_CGAL ${LIBIGL_TOPLEVEL_PROJECT})

# ディレクトリパスの設定 
# CMAKE_BINARY_DIRは，CMakeがビルドを行うために使用するディレクトリのパスを示す
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# CMake ヘルパー関数をインクルード
# これらのヘルパー関数は，libiglのビルドプロセスを簡素化するために使用される
include(igl_add_library)
include(igl_add_executable)
include(igl_copy_dll)
include(igl_include)
include(igl_install)
include(igl_target_sources)

# libiglのモジュールをCMakeのターゲットとして定義する
message(STATUS "=== include libigl modules ===")
igl_include(core)
igl_include(glfw)
igl_include(copyleft cgal)

# srcディレクトリを追加
add_subdirectory(src)

