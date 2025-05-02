#include "mesh.h"

// Constructor
migl::mesh::mesh( const int dim, const int num_vertices, const int num_faces )
    : m_dim(dim), m_vertex_matrix(Eigen::MatrixXd::Zero(num_vertices, dim)), m_face_matrix(Eigen::MatrixXi::Zero(num_faces, dim))
{
}

// Copy constructor
migl::mesh::mesh( const mesh& m )
    : m_dim(m.m_dim), m_vertex_matrix(m.m_vertex_matrix), m_face_matrix(m.m_face_matrix)
{
}

// Assignment operator
migl::mesh& migl::mesh::operator=( const mesh& m )
{
    if (this != &m)
    {
        m_dim = m.m_dim;
        m_vertex_matrix = m.m_vertex_matrix;
        m_face_matrix = m.m_face_matrix;
    }
    return *this;
}

// Output function
void migl::mesh::output() const
{
    std::cout << "Mesh Dimension: " << m_dim << std::endl;
    std::cout << "Vertex Matrix: " << std::endl << m_vertex_matrix << std::endl;
    std::cout << "Face Matrix: " << std::endl << m_face_matrix << std::endl;
}
