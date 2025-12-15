#ifndef MIGL_MESH_H
#define MIGL_MESH_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <iomanip>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/per_face_normals.h>
#include <Eigen/Dense>

namespace gmsh
{
    class geo;
}

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
            int num_vertices() const { return m_vertex_matrix.rows(); }
            int num_faces() const { return m_face_matrix.rows(); }
            const Eigen::MatrixXd& vertex_matrix() const { return m_vertex_matrix; }
            const Eigen::MatrixXd& V() const { return m_vertex_matrix; }
            const Eigen::RowVector3d V( const int i ) const { return m_vertex_matrix.row(i); }
            const double V( const int i, const int j ) const { return m_vertex_matrix(i,j); }
            const Eigen::MatrixXi& face_matrix() const { return m_face_matrix; }
            const Eigen::MatrixXi& F() const { return m_face_matrix; }
            const int F( const int i, const int j ) const { return m_face_matrix(i,j); }
            Eigen::Vector3d max_vertex() const;
            Eigen::Vector3d min_vertex() const;
            // Setters
            int& dim() { return m_dim; }
            Eigen::MatrixXd& vertex_matrix() { return m_vertex_matrix; }
            Eigen::MatrixXd& V() { return m_vertex_matrix; }
            //Eigen::RowVector3d& V( const int i ) { return m_vertex_matrix.row(i); }
            double& V( const int i, const int j ) { return m_vertex_matrix(i,j); }
            Eigen::MatrixXi& face_matrix() { return m_face_matrix; }
            Eigen::MatrixXi& F() { return m_face_matrix; }
            int& F( const int i, const int j ) { return m_face_matrix(i,j); }
            // Print function
            void output() const;
            std::ostream& output_vertex_matrix( std::ostream &out=std::cout, const int width=20, const int dec=10 ) const;
            std::ostream& output_face_matrix( std::ostream &out=std::cout, const int width=7 ) const;

            void orient_faces_consistently();
            bool is_outward_facing() const;
            void orient_faces_outward();
            void shift_vertices( const Eigen::Vector3d &shift );

            static void orient_faces_consistently( const Eigen::MatrixXi& F, Eigen::MatrixXi& F_out );
            static bool is_outward_facing( const Eigen::MatrixXd& V, const Eigen::MatrixXi& F );
            static mesh integrate_mesh( const std::vector<mesh>& meshes, std::vector<std::vector<int>>& face_indices_in_volumes );
            static gmsh::geo get_gmsh_geo( const std::vector<mesh>& meshes );
        
        private:
            static std::map<int,int> create_integrate_vertex_map( const mesh& base_map, const mesh& add_map );

        private:
            int m_dim;
            Eigen::MatrixXd m_vertex_matrix;
            Eigen::MatrixXi m_face_matrix;
    }; //mesh

} // namespace migl

#endif // MIGL_MESH_H