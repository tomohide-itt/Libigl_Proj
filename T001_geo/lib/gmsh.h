#ifndef GMESH_H
#define GMESH_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace gmsh
{
    class geo
    {
        public:
            static void read( const std::string& filename );

    };
}

#endif // GMESH_H