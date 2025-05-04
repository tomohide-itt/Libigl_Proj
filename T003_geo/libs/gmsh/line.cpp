#include "gmsh.h"

// Constructor for line class
gmsh::line::line(int tag, int start_point_tag, int end_point_tag)
    : m_tag(tag), m_start_point_tag(start_point_tag), m_end_point_tag(end_point_tag) {}

// Copy constructor for line class
gmsh::line::line(const line &l)
    : m_tag(l.m_tag), m_start_point_tag(l.m_start_point_tag), m_end_point_tag(l.m_end_point_tag) {}

// Assignment operator for line class
gmsh::line &gmsh::line::operator=(const line &l)
{
    if (this != &l)
    {
        m_tag = l.m_tag;
        m_start_point_tag = l.m_start_point_tag;
        m_end_point_tag = l.m_end_point_tag;
    }
    return *this;
}

// Output function for line class
void gmsh::line::output() const
{
    std::cout << "Line: " << m_tag << " (" << m_start_point_tag << ", " << m_end_point_tag << ")" << std::endl;
}