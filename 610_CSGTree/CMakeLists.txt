# 現在のディレクトリがトップレベルプロジェクトかどうかを確認
# このCMakeLists.txtがプロジェクトのトップレベルにある場合，LIBGL_PARENT_DIRは空になる
# そうでない場合，LIBGL_PARENT_DIRは親ディレクトリのパスになる
get_directory_property(LIBIGL_PARENT_DIR  PARENT_DIRECTORY)
# LIBIGL_PARENT_DIRが空であれば，LIBIGL_TOPLEVEL_PROJECTをONに設定
# そうでない場合，LIBIGL_TOPLEVEL_PROJECTをOFFに設定
if(NOT LIBIGL_PARENT_DIR)
    set(LIBIGL_TOPLEVEL_PROJECT ON)
else()
    set(LIBIGL_TOPLEVEL_PROJECT OFF)
endif()
message(STATUS "LIBIGL_TOPLEVEL_PROJECT: ${LIBIGL_TOPLEVEL_PROJECT}")


# CMakeのバージョンを確認
set(REQUIRED_CMAKE_VERSION "3.16.0")
# トップレベルプロジェクトであれば，CMakeの最低必要なバージョンを設定
# そうでない場合，CMakeのバージョンがREQUIRED_CMAKE_VERSION以下であれば，エラーメッセージを表示して終了
if(LIBIGL_TOPLEVEL_PROJECT)
    cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
else()
    if(${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
        message(FATAL_ERROR "CMake required version to build Libigl is ${REQUIRED_CMAKE_VERSION}")
    endif()
endif()


# ユーザが提供するデフォルトオプションファイル(LibiglOptions.cmake)がこのディレクトリに存在する場合，
# それをインクルードしてオプションを設定する
# これにより，CMakeLists.txtの他の部分で使用するオプションを定義できる
# `project()`の前であればインクルードすることで，C/C++コンパイラの設定（CMAKE_C_COMPILER, CMAKE_CXX_COMPILERなど）が有効になるため
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/LibiglOptions.cmake)
    message(STATUS "Using local options file: ${CMAKE_CURRENT_SOURCE_DIR}/LibiglOptions.cmake")
    include(${CMAKE_CURRENT_SOURCE_DIR}/LibiglOptions.cmake)
else()
    message( STATUS "No local options file found: ${CMAKE_CURRENT_SOURCE_DIR}/LibiglOptions.cmake")
endif()



# CMakeのオプションHUNTER_ENABLEDがONに設定されている場合，Hunterパッケージマネージャをセットアップする
# HUNTER_ENABLEDがOFFの場合，Hunterは使用されない（デフォルトではOFF）
# -DHUNTER_ENABLED=ONでCMakeを実行すると，Hunterが有効になる
# HunterはC++ライブラリの依存関係を管理するためのツールで，CMakeと統合されている
# Hunterを使用すると，CMakeLists.txtに依存関係を追加するだけで，Hunterが自動的に必要なライブラリをダウンロードしてビルドしてくれる
option(HUNTER_ENABLED "Enable Hunter package manager support" OFF)
if(HUNTER_ENABLED)
    include("cmake/misc/HunterGate.cmake")
    HunterGate(
        URL "https://github.com/cpp-pm/hunter/archive/v0.24.8.tar.gz"
        SHA1 "ca7838dded9a1811b04ffd56175f629e0af82d3d"
    )
endif()

################################################################################
# プロジェクトの名前とバージョンを設定
# プロジェクト名は610_CSGTree
# プロジェクトのバージョンは2.5.0
# プロジェクトのバージョンはlibiglのバージョンと同じにしている
project(610_CSGTree VERSION 2.5.0)


# CMakeのモジュール検索パス（CMAKE_MODULE_PATH）に複数のディレクトリを先頭に追加（PREPEND）
# これにより，CMakeは指定されたディレクトリ内のモジュールを優先的に検索する
# CMAKE_CURRENT_LIST_DIRは現在のCMakeLists.txtファイルのディレクトリの絶対パスを示す
list(PREPEND CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_LIST_DIR}/../libigl/cmake/
    ${CMAKE_CURRENT_LIST_DIR}/../libigl/cmake/igl
    ${CMAKE_CURRENT_LIST_DIR}/../libigl/cmake/find
    ${CMAKE_CURRENT_LIST_DIR}/../libigl/cmake/recipes/external
)
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")


# HUNTER_ENABLEDがONの場合，cmake/recipes/hunterディレクトリをCMAKE_MODULE_PATHに追加
# そうでない場合，cmake/recipes/externalディレクトリをCMAKE_MODULE_PATHに追加（この場合，同じパスが重複するが動作上は問題ない）
if(HUNTER_ENABLED)
    list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/../libigl/cmake/recipes/hunter)
else()
    list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/../libigl/cmake/recipes/external)
endif()

# CMakeのモジュールパスをグローバルプロパティ__igl_module_pathに設定
# これにより，libiglのモジュールを他のCMakeLists.txtファイルからも参照できるようになる
# __igl_module_pathはlibiglのモジュールを検索するためのパスを格納するプロパティ
set_property(GLOBAL PROPERTY __igl_module_path ${CMAKE_MODULE_PATH})


# トップレベルなら，LIBIGL_DEFAULT_CGALをONに設定
set(LIBIGL_DEFAULT_CGAL ${LIBIGL_TOPLEVEL_PROJECT})
# MATLABに対する設定
# 特定バージョンの MATLAB に対して，対応する内部バージョン番号（APIバージョン）をマッピング
# これは MATLAB サポートを検出するために使われることがあり，CMake の FindMatlab.cmake などが
# 対応していないバージョンでも見つけられるようにするため
set(MATLAB_ADDITIONAL_VERSIONS
    "R2023b=10.4"
    "R2023a=10.4"
    "R2022b=10.3"
    "R2022a=10.2"
    "R2021b=10.1"
    "R2021a=10.0"
)
# トップレベルなら，LIBIGL_DEFAULT_MATLABをONに設定
set(LIBIGL_DEFAULT_MATLAB ${LIBIGL_TOPLEVEL_PROJECT})
# トップレベルなら，LIBIGL_DEFAULT_MOSEKをONに設定
set(LIBIGL_DEFAULT_MOSEK ${LIBIGL_TOPLEVEL_PROJECT})
# トップレベルプロジェクトとしてビルドしている場合に限り，
# MATLAB と MOSEK の自動検出を試みて，見つからなければそれぞれのモジュールを無効化する処理
if(LIBIGL_TOPLEVEL_PROJECT)
    # MATLABの検出
    # QUIETにより，見つからなくてもエラーメッセージを表示しない
    # find_packageは，自身のプロジェクトに属していないライブラリを自動検索するコマンド
    # モジュールモード，コンフィグモードの順で探す．
    # モジュールモードでは，Find<Package>.cmakeを探し，コンフィグモードでは，<Package>Config.cmakeを探す
    # FindMatlab.cmake を使って、MATLAB を探す
    # MatlabConfig.cmake 
    find_package(Matlab QUIET)
    # MATLABが見つからない場合，LIBIGL_DEFAULT_MATLABをOFFに設定して警告を表示
    if(NOT Matlab_FOUND)
        set(LIBIGL_DEFAULT_MATLAB OFF)
        message(WARNING "Matlab not found, disabling igl_restricted::matlab module.")
    endif()
    # MOSEK（数理最適化ソルバ）の検出
    # FindMOESK.cmake を使って、MOSEK を探す
    find_package(MOSEK QUIET)
    # MOSEKが見つからない場合，LIBIGL_DEFAULT_MOSEKをOFFに設定して警告を表示
    if(NOT MOSEK_FOUND)
        set(LIBIGL_DEFAULT_MOSEK OFF)
        message(WARNING "Mosek not found, disabling igl_restricted::mosek module.")
    endif()
endif()

# トップレベルなら，LIBIGL_DEFAULT_COMISOをONに設定
# ただし，APPLEプラットフォームでない場合やBLASライブラリが見つからなければOFFに設定
set(LIBIGL_DEFAULT_COMISO ${LIBIGL_TOPLEVEL_PROJECT})
if(LIBIGL_TOPLEVEL_PROJECT AND (NOT APPLE) AND UNIX)
  find_package(BLAS QUIET)
    if(NOT BLAS_FOUND)
        set(LIBIGL_DEFAULT_COMISO OFF)
        message(WARNING "BLAS not found, disabling igl_copyleft::comiso module.")
    endif()
endif()

# トップレベルなら，LIBIGL_BUILD_TESTS, LIBIGL_BUILD_TUTORIALS, LIBIGL_INSTALLをONに設定
option(LIBIGL_BUILD_TESTS      "Build libigl unit test"                ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_BUILD_TUTORIALS  "Build libigl tutorial"                 ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_INSTALL          "Enable installation of libigl targets" ${LIBIGL_TOPLEVEL_PROJECT})


# LIBIGL_USE_STATIC_LIBRARYがONの場合，libiglを静的ライブラリとしてビルドする
# ただし，LIBIGL_TOPLEVEL_PROJECTがONの場合のみ
# これは，libiglがヘッダオンリーライブラリであるため，デフォルトではOFFに設定されている
# ONにすることで，複数のバイナリの生成が高速化される
# ただし，初期のコンパイル時間は長くなる
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" ${LIBIGL_TOPLEVEL_PROJECT})


# これらのモジュールは，MPL2ライセンスで利用可能であり，その依存関係は許可されたライセンスまたはパブリックドメインライセンスで利用可能
# これらのモジュールは，libiglのコア機能を拡張するために使用される
option(LIBIGL_EMBREE           "Build target igl::embree"           ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_GLFW             "Build target igl::glfw"             ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_IMGUI            "Build target igl::imgui"            ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_OPENGL           "Build target igl::opengl"           ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_STB              "Build target igl::stb"              ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_PREDICATES       "Build target igl::predicates"       ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_SPECTRA          "Build target igl::spectra"          ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_XML              "Build target igl::xml"              ${LIBIGL_TOPLEVEL_PROJECT})


# これらのモジュールは，GPLライセンスで利用可能であり，その依存関係はコピーレフトライセンスで利用可能
option(LIBIGL_COPYLEFT_CORE    "Build target igl_copyleft::core"    ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_COPYLEFT_CGAL    "Build target igl_copyleft::cgal"    ${LIBIGL_DEFAULT_CGAL})
option(LIBIGL_COPYLEFT_COMISO  "Build target igl_copyleft::comiso"  ${LIBIGL_TOPLEVEL_PROJECT})
option(LIBIGL_COPYLEFT_TETGEN  "Build target igl_copyleft::tetgen"  ${LIBIGL_TOPLEVEL_PROJECT})


# これらのモジュールは，MPL2ライセンスで利用可能であり，その依存関係は非商用ライセンスまたはプロプライエタリライセンスで利用可能
option(LIBIGL_RESTRICTED_MATLAB   "Build target igl_restricted::matlab"   ${LIBIGL_DEFAULT_MATLAB})
option(LIBIGL_RESTRICTED_MOSEK    "Build target igl_restricted::mosek"    ${LIBIGL_DEFAULT_MOSEK})
option(LIBIGL_RESTRICTED_TRIANGLE "Build target igl_restricted::triangle" ${LIBIGL_TOPLEVEL_PROJECT})


# GLFWはヘッドレスCIマシンでは動作しないため，テストを実行しない（またはビルドしない）
# 不幸なことに，ヘッドレスmacマシンではglfwがクラッシュするのではなくハングアップするようで，
# 実行時にキャッチするのが難しい
option(LIBIGL_GLFW_TESTS       "Build igl::glfw tests"             ${LIBIGL_TOPLEVEL_PROJECT})

# LIBIGL_WARNINGS_AS_ERRORSをOFFに設定
option(LIBIGL_WARNINGS_AS_ERRORS "Turn on many warnings and treat as errors" OFF)

# ディレクトリパスの設定 
# CMAKE_BINARY_DIRは，CMakeがビルドを行うために使用するディレクトリのパスを示す
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")


# CMake ヘルパー関数をインクルード
# cmake/igl_add_library.cmakeは，ライブラリを追加するためのヘルパー関数を定義している
# cmake/igl_add_test.cmakeは，テストを追加するためのヘルパー関数を定義している
# cmake/igl_add_tutorial.cmakeは，チュートリアルを追加するためのヘルパー関数を定義している
# cmake/igl_copy_dll.cmakeは，DLLをコピーするためのヘルパー関数を定義している
# cmake/igl_include.cmakeは，ヘッダーファイルをインクルードするためのヘルパー関数を定義している
# cmake/igl_install.cmakeは，インストールするためのヘルパー関数を定義している
# cmake/igl_target_sources.cmakeは，ターゲットのソースファイルを追加するためのヘルパー関数を定義している
# これらのヘルパー関数は，libiglのビルドプロセスを簡素化するために使用される
include(igl_add_library)
include(igl_add_test)
include(igl_add_tutorial)
include(igl_copy_dll)
include(igl_include)
include(igl_install)
include(igl_target_sources)


# ルートレベルでユニットテストを有効にする
# ユニットテストは，個々のモジュールや関数が正しく動作するかどうかを確認するためのテスト
if(LIBIGL_BUILD_TESTS)
    include(CTest)
endif()


# これにより，libiglのモジュールをCMakeのターゲットとして定義する
include(libigl)

#[[
# tutorialディレクトリを追加
# これは，libiglのチュートリアルをビルドするためのディレクトリ
# チュートリアルは，libiglの使い方を学ぶためのサンプルコードやドキュメントを含む
if(LIBIGL_BUILD_TUTORIALS)
    add_subdirectory(tutorial)
endif()

################################################################################
# Install CMake config files
################################################################################

if(LIBIGL_INSTALL)
    include(GNUInstallDirs)
    set(project_config_in "${PROJECT_SOURCE_DIR}/cmake/igl/libigl-config.cmake.in")
    set(project_config_out "${CMAKE_CURRENT_BINARY_DIR}/libigl-config.cmake")
    set(config_targets_base "LibiglConfigTargets")
    set(version_config_file "${CMAKE_CURRENT_BINARY_DIR}/LibiglConfigVersion.cmake")
    set(export_dest_dir "${CMAKE_INSTALL_LIBDIR}/cmake/igl")

    foreach(suffix IN ITEMS "") #"_restricted" "_copyleft")
        install(EXPORT LibiglTargets${suffix}
            DESTINATION ${export_dest_dir}
            NAMESPACE igl${suffix}::
            FILE ${config_targets_base}${suffix}.cmake
            COMPONENT LibiglDevelopment
        )
    endforeach()

    include(CMakePackageConfigHelpers)
    configure_package_config_file(
            "${project_config_in}"
            "${project_config_out}"
        INSTALL_DESTINATION
            ${CMAKE_INSTALL_DATAROOTDIR}/libigl/cmake
    )

    write_basic_package_version_file("${version_config_file}" COMPATIBILITY SameMajorVersion)
    install(FILES "${project_config_out}" "${version_config_file}" DESTINATION "${export_dest_dir}")
endif()
]]
