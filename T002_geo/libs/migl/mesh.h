#ifndef MIGL_MESH_H
#define MIGL_MESH_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>

namespace migl
{
    class mesh
    {
        public:
            mesh() : m_dim(0) {}
            mesh( const int dim, const int num_vertices, const int num_faces );
            mesh( const mesh& m );
            mesh& operator=( const mesh& m );
            // Getters
            int dim() const { return m_dim; }
            const Eigen::MatrixXd& vertex_matrix() const { return m_vertex_matrix; }
            const Eigen::MatrixXi& face_matrix() const { return m_face_matrix; }
            // Setters
            int& dim() { return m_dim; }
            Eigen::MatrixXd& vertex_matrix() { return m_vertex_matrix; }
            Eigen::MatrixXi& face_matrix() { return m_face_matrix; }
            // Print function
            void output() const;

        private:
            int m_dim;
            Eigen::MatrixXd m_vertex_matrix;
            Eigen::MatrixXi m_face_matrix;
    }; //mesh

} // namespace migl

#endif // MIGL_MESH_H