set( LIB_NAME "gmsh" )

add_library( ${LIB_NAME} STATIC
    gmsh.h
    geo.cpp
    point.cpp
    line.cpp
    surface.cpp
    volume.cpp
)