set( LIB_NAME "migl" )

set( HDRS
    mesh.h
)

set( SRCS 
    mesh.cpp
)

add_library( ${LIB_NAME} STATIC ${HDRS} ${SRCS} )

target_compile_features(${LIB_NAME} PUBLIC cxx_std_17)

target_link_libraries(${LIB_NAME} PUBLIC
    igl::core
)

target_include_directories( ${LIB_NAME} PRIVATE
    ${LIB_DIR}/gmsh
)

target_compile_definitions( ${LIB_NAME} PUBLIC DEBUG_PATH=\"${CMAKE_BINARY_DIR}/debug\" )