#include "gmsh.h"

// Constructor for surface class
gmsh::surface::surface(int tag, const std::vector<int>& line_tags )
    : m_tag(tag), m_line_tags(line_tags) {}

// Copy constructor for surface class
gmsh::surface::surface(const surface &s)
    : m_tag(s.m_tag), m_line_tags(s.m_line_tags) {}

// Assignment operator for surface class
gmsh::surface &gmsh::surface::operator=(const surface &s)
{
    if (this != &s)
    {
        m_tag = s.m_tag;
        m_line_tags = s.m_line_tags;
    }
    return *this;
}

// Output function for surface class
void gmsh::surface::output() const
{
    std::cout << "Surface: " << m_tag << " (";
    for (size_t i = 0; i < m_line_tags.size(); ++i)
    {
        std::cout << m_line_tags[i];
        if (i < m_line_tags.size() - 1)
        {
            std::cout << ", ";
        }
    }
    std::cout << ")" << std::endl;
}