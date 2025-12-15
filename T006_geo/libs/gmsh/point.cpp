#include "gmsh.h"

// Constructor for point class
gmsh::point::point(double x, double y, double z, int tag)
    : m_x(x), m_y(y), m_z(z), m_tag(tag) {}

// Copy constructor for point class
gmsh::point::point(const point &p)
    : m_x(p.m_x), m_y(p.m_y), m_z(p.m_z), m_tag(p.m_tag) {}

// Assignment operator for point class
gmsh::point &gmsh::point::operator=(const point &p)
{
    if (this != &p)
    {
        m_x = p.m_x;
        m_y = p.m_y;
        m_z = p.m_z;
        m_tag = p.m_tag;
    }
    return *this;
}

// Output function for point class
void gmsh::point::output() const
{
    std::cout << "Point: " << m_tag << " (" << m_x << ", " << m_y << ", " << m_z << ")" << std::endl;
}