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
    std::cout << "Vertex Matrix: " << std::endl << m_vertex_matrix << std::endl;
    std::cout << "Face Matrix: " << std::endl << m_face_matrix << std::endl;
}

// Orient faces consistently (member function)
void migl::mesh::orient_faces_consistently()
{
    Eigen::MatrixXi F_out;
    orient_faces_consistently(m_face_matrix, F_out);
    m_face_matrix = F_out;
}

// Check if the faces are outward facing (member function)
bool migl::mesh::is_outward_facing() const
{
    return is_outward_facing(m_vertex_matrix, m_face_matrix);
}

// Orient faces outward (member function)
void migl::mesh::orient_faces_outward()
{
    this->orient_faces_consistently();
    if( !this->is_outward_facing() )
    {
        Eigen::MatrixXi F_out = this->face_matrix();
        for( int i=0; i<F_out.rows(); i++ )
        {
            std::swap( F_out(i,0), F_out(i,1) );
        }
        this->face_matrix() = F_out;
    }
}

// Orient faces consistently (static function)
void migl::mesh::orient_faces_consistently( const Eigen::MatrixXi& F, Eigen::MatrixXi& F_out )
{
    F_out = F;
    
    Eigen::MatrixXi TT;
    igl::triangle_triangle_adjacency(F, TT);
    //+++
    //std::cout << "TT: " << std::endl << TT << std::endl;
    //---
    int num_faces = F.rows();
    std::vector<bool> visited(num_faces, false);
    std::queue<int> q;

    q.push(0);
    visited[0] = true;

    while (!q.empty())
    {
        int fid = q.front();
        //+++
        //std::cout << "fid: " << fid << std::endl;
        //---
        q.pop();

        for (int i = 0; i < 3; ++i)
        {
            int nbr = TT(fid, i);
            if (nbr < 0 || visited[nbr]) continue;

            int v0 = F_out(fid,i);
            int v1 = F_out(fid,(i+1)%3);

            int match_count = 0;
            int match_dir = 0; //1:same, -1:reverse

            for( int j=0; j<3; j++ )
            {
                int nv0 = F_out(nbr,j);
                int nv1 = F_out(nbr,(j+1)%3);
                if( v0 == nv0 && v1 == nv1 )
                {
                    match_count++;
                    match_dir = 1;
                }
                else if( v0 == nv1 && v1 == nv0 )
                {
                    match_count++;
                    match_dir = -1;
                }
            }

            if( match_count == 1 )
            {
                if( match_dir == 1 )
                {
                    std::swap( F_out(nbr,0), F_out(nbr,1) );
                    std::swap( TT(nbr,1), TT(nbr,2));
                }
            }
            
            visited[nbr] = true;
            q.push(nbr);
            //+++
            //std::cout << "  nbr: " << nbr << " match_count: " << match_count << " match_dir: " << match_dir << std::endl;
            //std::cout << "  F_out.row(nbr): " << F_out.row(nbr) << std::endl;
            //---
        }
        //+++
        //std::cout << "  q: ";
        //std::queue<int> q_copy = q;
        //while (!q_copy.empty())
        //{
        //    std::cout << q_copy.front() << " ";
        //    q_copy.pop();
        //}
        //std::cout << std::endl;
        //---
    }
}

// Check if the faces are outward facing
bool migl::mesh::is_outward_facing( const Eigen::MatrixXd& V, const Eigen::MatrixXi& F )
{
    // Center of the bounding box
    Eigen::RowVector3d center = V.colwise().mean();

    // Compute the normals
    Eigen::MatrixXd N;
    igl::per_face_normals(V, F, N);

    int outward_count = 0;
    int inward_count = 0;

    for (int i = 0; i < F.rows(); ++i)
    {
        Eigen::RowVector3d p0 = V.row(F(i, 0));
        Eigen::RowVector3d p1 = V.row(F(i, 1));
        Eigen::RowVector3d p2 = V.row(F(i, 2));

        // Compute the centroid of the face
        Eigen::RowVector3d centroid = (p0 + p1 + p2) / 3.0;
        // Compute the direction from the centroid to the center
        Eigen::RowVector3d dir = center - centroid;
        // Check if the dot product is positive or negative
        double dot_product = dir.dot(N.row(i));
        if (dot_product < 0)
        {
            outward_count++;
        }
        else
        {
            inward_count++;
        }
    }
    return (outward_count >= inward_count);
}