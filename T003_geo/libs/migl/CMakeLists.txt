set( LIB_NAME "migl" )

set( HDRS
    mesh.h
)

set( SRCS 
    mesh.cpp
)

add_library( ${LIB_NAME} STATIC ${HDRS} ${SRCS} )

target_link_libraries(${LIB_NAME} PUBLIC
    igl::core
)