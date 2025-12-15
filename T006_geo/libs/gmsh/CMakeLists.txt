set( LIB_NAME "gmsh" )

set( HDRS
    gmsh.h
)

set( SRCS 
    geo.cpp
    point.cpp
    line.cpp
    surface.cpp
    volume.cpp
)

add_library( ${LIB_NAME} STATIC ${HDRS} ${SRCS} )

target_compile_features(${LIB_NAME} PUBLIC cxx_std_17)

target_link_libraries(${LIB_NAME} PUBLIC
    igl::core
    migl
)

target_include_directories( ${LIB_NAME} PRIVATE
    ${LIB_DIR}/migl
)

message( STATUS "LIB_DIR: ${LIB_DIR}" )