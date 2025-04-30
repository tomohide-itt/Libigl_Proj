#必要なCMakeのバージョンを指定
cmake_minimum_required(VERSION 3.16)

####################################################################################
# プロジェクトの名前とバージョンを設定
project(610_CSGTree VERSION 2.5.0)

# CMakeのモジュール検索パス（CMAKE_MODULE_PATH）に複数のディレクトリを先頭に追加（PREPEND）
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libiglのモジュールをCMakeのターゲットとして定義する
message(STATUS "=== include libigl modules ===")
include(libigl)

# モジュールを有効化する
igl_include(core)
igl_include(glfw)
igl_include(copyleft cgal)

# ソースファイルを集める 
file(GLOB SRC_FILES *.cpp)

# 実行ファイルを生成する 
add_executable(${PROJECT_NAME} ${SRC_FILES})

# 必要なライブラリをリンクする 
target_link_libraries(${PROJECT_NAME} PUBLIC 
  igl::core
)

# データへのパスを指定する 
target_compile_definitions( ${PROJECT_NAME} PRIVATE TUTORIAL_SHARED_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data\" )