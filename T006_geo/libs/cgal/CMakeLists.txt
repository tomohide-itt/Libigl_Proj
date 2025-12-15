set( LIB_NAME "cgal" )

set( HDRS
    ${LIB_NAME}.h
)

set( SRCS 
    ${LIB_NAME}.cpp
)

add_library( ${LIB_NAME} STATIC ${HDRS} ${SRCS} )

target_compile_features(${LIB_NAME} PUBLIC cxx_std_17)

target_link_libraries(${LIB_NAME} PUBLIC
    igl_copyleft::cgal
)
#
#target_include_directories( ${LIB_NAME} PRIVATE
#    ${LIB_DIR}/gmsh
#)