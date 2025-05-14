#include "gmsh.h"

// Constructor for volume class
gmsh::volume::volume(int tag, const std::vector<int>& surface_tags )
    : m_tag(tag), m_surface_tags(surface_tags) {}

// Copy constructor for volume class
gmsh::volume::volume(const volume &v)
    : m_tag(v.m_tag), m_surface_tags(v.m_surface_tags) {}

// Assignment operator for volume class
gmsh::volume &gmsh::volume::operator=(const volume &v)
{
    if (this != &v)
    {
        m_tag = v.m_tag;
        m_surface_tags = v.m_surface_tags;
    }
    return *this;
}

// Output function for volume class 
void gmsh::volume::output() const
{
    std::cout << "Volume: " << m_tag << " (";
    for (size_t i = 0; i < m_surface_tags.size(); ++i)
    {
        std::cout << m_surface_tags[i];
        if (i < m_surface_tags.size() - 1)
        {
            std::cout << ", ";
        }
    }
    std::cout << ")" << std::endl;
}