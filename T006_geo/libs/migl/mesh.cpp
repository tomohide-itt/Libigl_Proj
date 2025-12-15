#include "mesh.h"
#include "gmsh.h"

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

// Get the maximum vertex (member function)
Eigen::Vector3d migl::mesh::max_vertex() const
{
    Eigen::Vector3d max_vertex = m_vertex_matrix.colwise().maxCoeff();
    return max_vertex;
}

// Get the minimum vertex (member function)
Eigen::Vector3d migl::mesh::min_vertex() const
{
    Eigen::Vector3d min_vertex = m_vertex_matrix.colwise().minCoeff();
    return min_vertex;
}

// Output function
void migl::mesh::output() const
{
    std::cout << "Vertex Matrix: " << std::endl << m_vertex_matrix << std::endl;
    std::cout << "Face Matrix: " << std::endl << m_face_matrix << std::endl;
}

//
std::ostream& migl::mesh::output_vertex_matrix( std::ostream &out, const int width, const int dec ) const
{
    out << std::scientific << std::setprecision(dec);
    for( int i=0; i<this->m_vertex_matrix.rows(); i++ )
    {
        out << std::setw(10) << i;
        out << std::setw(width) << this->m_vertex_matrix(i,0);
        out << std::setw(width) << this->m_vertex_matrix(i,1);
        out << std::setw(width) << this->m_vertex_matrix(i,2);
        out << std::endl;
    }
    out << std::defaultfloat;
    return out;
}
//
std::ostream& migl::mesh::output_face_matrix( std::ostream &out, const int width ) const
{
    for( int i=0; i<this->m_face_matrix.rows(); i++ )
    {
        out << std::setw(10) << i;
        out << std::setw(width) << this->m_face_matrix(i,0);
        out << std::setw(width) << this->m_face_matrix(i,1);
        out << std::setw(width) << this->m_face_matrix(i,2);
        out << std::endl;
    }
    return out;
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

void migl::mesh::shift_vertices( const Eigen::Vector3d &shift )
{
    for( int i=0; i<m_vertex_matrix.rows(); i++ )
    {
        m_vertex_matrix(i,0) += shift(0);
        m_vertex_matrix(i,1) += shift(1);
        m_vertex_matrix(i,2) += shift(2);
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

// Integrate multiple meshes into one
migl::mesh migl::mesh::integrate_mesh( const std::vector<migl::mesh>& meshes, std::vector<std::vector<int>>& face_indices_in_volumes )
{
//+++
    bool debug = true;
    std::ofstream fdeb(DEBUG_PATH "/integrate_mesh.txt");
//---
    //初期の統合メッシュにmeshes[0]を代入する
    migl::mesh integrated_mesh = meshes[0];
//+++
    /*
    if( debug )
    {
        fdeb << "integrating mesh[" << 0 << "] ------------------" << std::endl;
        fdeb << "base_mesh V:" << std::endl << integrated_mesh.vertex_matrix() << std::endl;
        fdeb << "base_mesh F:" << std::endl << integrated_mesh.face_matrix() << std::endl;
    }
    */
//---
    //face_indices_in_volumes: 各Volumeを形成するfaceのindexを格納するための変数
    //その初期設定
    face_indices_in_volumes.resize(meshes.size());
    for( int i=0; i<meshes[0].face_matrix().rows(); i++ )
    {
        face_indices_in_volumes[0].push_back(i);
    }
    //
    for( int i=1; i<meshes.size(); i++ )
    {
        std::cout << "integrating mesh[" << i << "] / " << meshes.size()-1 << std::endl;
//+++
        debug = false;
        //if( i==7 ) debug = true;
        if( debug )
        {
            fdeb << "integrating mesh[" << i << "] ------------------" << std::endl;
            fdeb << "base_mesh V:" << std::endl;
            integrated_mesh.output_vertex_matrix(fdeb);
            fdeb << "base_mesh F:" << std::endl;
            integrated_mesh.output_face_matrix(fdeb);
            fdeb << "mesh[" << i << "] V:" << std::endl;
            meshes[i].output_vertex_matrix(fdeb);
            fdeb << "mesh[" << i << "] F:" << std::endl;
            meshes[i].output_face_matrix(fdeb);
        }
//---
        //メッシュ[i]の頂点マトリクスを統合するため，
        //メッシュ[i]の頂点マトリクスのindexをkey
        //統合後の頂点マトリクスのindexをvalueとするマップを作成する
        std::map<int,int> vertex_map = create_integrate_vertex_map( integrated_mesh, meshes[i] );
//+++
        if( debug )
        {
            fdeb << "vertex_map: " << std::endl;
            for( const auto &[key,val] : vertex_map )
            {
                fdeb << key << " -> " << val << std::endl;
            }
        }
//---
        //統合メッシュの頂点数
        int base_mesh_num_vertices = integrated_mesh.vertex_matrix().rows();
        //メッシュ[i]の頂点を統合するとき，新たに追加する頂点数(add_mesh_num_vertices)を数える
        int add_mesh_num_vertices = 0;
        for( const auto& [add_vindex, vindex] : vertex_map )
        {
            if( vindex >= base_mesh_num_vertices )
            {
                add_mesh_num_vertices++;
            }
        }
        //統合後の頂点数
        int num_vertices = base_mesh_num_vertices + add_mesh_num_vertices;
        //統合メッシュのV
        Eigen::MatrixXd base_V = integrated_mesh.vertex_matrix();
        //統合メッシュに追加するV
        Eigen::MatrixXd add_V(add_mesh_num_vertices, 3);
        //mesh[i]から統合メッシュのVにはない頂点を選択して，add_Vに登録する
        int index = 0;
        for( const auto& [add_vindex, vindex] : vertex_map )
        {
            if( vindex >= base_mesh_num_vertices )
            {
                add_V.row(index++) = meshes[i].vertex_matrix().row(add_vindex);
            }
        }
        //統合メッシュの頂点マトリクスのサイズを更新し，メッシュ[i]の頂点を追加
        integrated_mesh.vertex_matrix().resize(num_vertices, 3);
        integrated_mesh.vertex_matrix() << base_V, add_V;
        //
        //統合メッシュの面の数
        int base_mesh_num_faces = integrated_mesh.face_matrix().rows();
        int add_mesh_num_faces = 0;
        //
        //メッシュ[i]のj番目の面を形成する頂点のすべてが統合メッシュの頂点に含まれている場合，trueを
        //そうでなければfalseを返す　ラムダ関数
        auto out_of_range = [&](int j) ->bool
        {
//+++
            //if (debug)
            //{
            //    fdeb << "meshes[" << i << "]::face[" << j << "]:( ";
            //    fdeb << meshes[i].face_matrix()(j, 0) << ", ";
            //    fdeb << meshes[i].face_matrix()(j, 1) << ", ";
            //    fdeb << meshes[i].face_matrix()(j, 2) << " )" << std::endl;
            //}
//---
            bool ret = vertex_map.at(meshes[i].face_matrix()(j,0)) >= base_mesh_num_vertices ||
                   vertex_map.at(meshes[i].face_matrix()(j,1)) >= base_mesh_num_vertices ||
                   vertex_map.at(meshes[i].face_matrix()(j,2)) >= base_mesh_num_vertices;
//+++
            //if( debug )
            //{
            //    fdeb << "converted to vertex index after integration:( ";
            //    fdeb << vertex_map.at(meshes[i].face_matrix()(j,0)) << ", ";
            //    fdeb << vertex_map.at(meshes[i].face_matrix()(j,1)) << ", ";
            //    fdeb << vertex_map.at(meshes[i].face_matrix()(j,2)) << " ) >= " << base_mesh_num_vertices << " ?" << std::endl;
            //    fdeb << "out_of_range: " << std::boolalpha << ret << std::endl;
            //}
//---
            return ret;
        };
        //
        //統合メッシュのF
        Eigen::MatrixXi base_F = integrated_mesh.face_matrix();
        //メッシュ[i]の面を統合するとき，新たに追加する面の数(add_mesh_num_faces)を数える
        for( int j=0; j<meshes[i].face_matrix().rows(); j++ )
        {
            if( out_of_range(j) )
            {
                add_mesh_num_faces++;
            }
            else
            {
                bool found_same_face = false;
                for( int k=0; k<base_mesh_num_faces; k++ )
                {
                    auto match_vertex_index = [&]( int ind0, int ind1, int ind2 )->bool
                    {
                        return base_F(k,0) == vertex_map.at(meshes[i].face_matrix()(j,ind0)) &&
                               base_F(k,1) == vertex_map.at(meshes[i].face_matrix()(j,ind1)) &&
                               base_F(k,2) == vertex_map.at(meshes[i].face_matrix()(j,ind2));
                    };
                    if( match_vertex_index(0,1,2) || match_vertex_index(0,2,1) ||
                        match_vertex_index(1,0,2) || match_vertex_index(1,2,0) ||
                        match_vertex_index(2,0,1) || match_vertex_index(2,1,0) )
                    {
                        found_same_face = true;
                        break;
                    }
                }
                //３つの頂点はすでに統合メッシュに存在しているが，面はまだ存在していない場合
                if( !found_same_face )
                {
                    add_mesh_num_faces++;
                }
            }
        }
        //統合後の面の数
        int num_faces = base_mesh_num_faces + add_mesh_num_faces;
        //統合メッシュに追加するF
        Eigen::MatrixXi add_F(add_mesh_num_faces, 3);
        //mesh[i]から統合メッシュのFにはない面を選択して，add_Fに登録する
        index = 0;
        for( int j=0; j<meshes[i].face_matrix().rows(); j++ )
        {
            if( out_of_range(j) )
            {
                face_indices_in_volumes[i].push_back( index + base_mesh_num_faces );
                add_F.row(index++) << vertex_map.at(meshes[i].face_matrix()(j,0)),
                                      vertex_map.at(meshes[i].face_matrix()(j,1)),
                                      vertex_map.at(meshes[i].face_matrix()(j,2));
            }
            else
            {
                bool found_same_face = false;
                for( int k=0; k<base_mesh_num_faces; k++ )
                {
                    auto match_vertex_index = [&]( int ind0, int ind1, int ind2 )->bool
                    {
                        return base_F(k,0) == vertex_map.at(meshes[i].face_matrix()(j,ind0)) &&
                               base_F(k,1) == vertex_map.at(meshes[i].face_matrix()(j,ind1)) &&
                               base_F(k,2) == vertex_map.at(meshes[i].face_matrix()(j,ind2));
                    };
                    if( match_vertex_index(0,1,2) || match_vertex_index(0,2,1) ||
                        match_vertex_index(1,0,2) || match_vertex_index(1,2,0) ||
                        match_vertex_index(2,0,1) || match_vertex_index(2,1,0) )
                    {
                        face_indices_in_volumes[i].push_back(k);
                        found_same_face = true;
                        break;
                    }
                }
                //３つの頂点はすでに統合メッシュに存在しているが，面はまだ存在していない場合
                if( !found_same_face )
                {
                    face_indices_in_volumes[i].push_back( index + base_mesh_num_faces );
                    add_F.row(index++) << vertex_map.at(meshes[i].face_matrix()(j,0)),
                                          vertex_map.at(meshes[i].face_matrix()(j,1)),
                                          vertex_map.at(meshes[i].face_matrix()(j,2));
                }
            }
        }
        integrated_mesh.face_matrix().resize(num_faces, 3);
        integrated_mesh.face_matrix() << base_F, add_F;
    }
    //+++
    /*
    std::cout << "integrated_mesh.vertex_matrix(): " << std::endl;
    std::cout << integrated_mesh.vertex_matrix() << std::endl;
    std::cout << "integrated_mesh.face_matrix(): " << std::endl;
    std::cout << integrated_mesh.face_matrix() << std::endl;
    std::cout << "indices_in_volumes: " << std::endl;
    for( int i=0; i<face_indices_in_volumes.size(); i++ )
    {
        std::cout << "  face_indices_in_volumes[" << i << "]: ";
        for( int j=0; j<face_indices_in_volumes[i].size(); j++ )
        {
            std::cout << face_indices_in_volumes[i][j] << " ";
        }
        std::cout << std::endl;
    }
    */
    //---
//+++
    fdeb.close();
//---
    return integrated_mesh;
}

// Create the vertex map for integration
std::map<int,int> migl::mesh::create_integrate_vertex_map( const mesh& base_map, const mesh& add_map )
{
    std::map<int,int> vertex_map;
    Eigen::Vector3d base_max = base_map.max_vertex();
    Eigen::Vector3d base_min = base_map.min_vertex();
    Eigen::Vector3d add_max = add_map.max_vertex();
    Eigen::Vector3d add_min = add_map.min_vertex();
    //
    int num_base_vertices = base_map.vertex_matrix().rows();
    int num_add_vertices = add_map.vertex_matrix().rows();
    int add_index = num_base_vertices;
    //
    for( int i=0; i<num_add_vertices; i++ )
    {
        Eigen::RowVector3d pt = add_map.vertex_matrix().row(i);
        if( pt(0) < base_min(0) || pt(0) > base_max(0) ||
            pt(1) < base_min(1) || pt(1) > base_max(1) ||
            pt(2) < base_min(2) || pt(2) > base_max(2) )
        {
            vertex_map[i] = add_index++;
        }
        else
        {
            bool found = false;
            for( int j=0; j<num_base_vertices; j++ )
            {
                Eigen::RowVector3d base_pt = base_map.vertex_matrix().row(j);
                if( pt.isApprox(base_pt, 1e-10) )
                {
                    vertex_map[i] = j;
                    found = true;
                    break;
                }
            }
            if( !found ) vertex_map[i] = add_index++;
        }
    }
    return vertex_map;
}

// Get the GMSH geo object from the meshes
gmsh::geo migl::mesh::get_gmsh_geo( const std::vector<migl::mesh>& meshes )
{
    gmsh::geo geo;
    std::vector<std::vector<int>> face_indices_in_volumes;
    migl::mesh mesh = integrate_mesh( meshes, face_indices_in_volumes );
    //
    std::vector<gmsh::point> points;
    for( int i=0; i<mesh.vertex_matrix().rows(); i++ )
    {
        gmsh::point pt;
        pt.tag() = i+1;
        pt.x() = mesh.vertex_matrix()(i,0);
        pt.y() = mesh.vertex_matrix()(i,1);
        pt.z() = mesh.vertex_matrix()(i,2);
        points.push_back(pt);
    }
    //
    std::vector<gmsh::line> lines;
    std::vector<gmsh::surface> surfaces;
    int new_ln_tag = 1;
    std::map<int,std::vector<int>> map_firstp2ln;
    for( int i=0; i<mesh.face_matrix().rows(); i++ )
    {
        gmsh::surface sf;
        sf.tag() = i+1;
        std::array<int,3> pt_tags = { mesh.face_matrix()(i,0)+1, mesh.face_matrix()(i,1)+1, mesh.face_matrix()(i,2)+1 };
        for( int j=0; j<3; j++ )
        {
            bool reverse = false;
            int start_pt_tag = pt_tags[j];
            int end_pt_tag = pt_tags[(j+1)%3];
            if( pt_tags[j] > pt_tags[(j+1)%3] )
            {
                std::swap(start_pt_tag, end_pt_tag);
                reverse = true;
            }
            std::vector<int> exist_ln_tags = map_firstp2ln[start_pt_tag];
            int ln_tag = -1;
            for( const auto &exist_ln_tag : exist_ln_tags )
            {
                if( lines[exist_ln_tag].end_point_tag() == end_pt_tag )
                {
                    ln_tag = exist_ln_tag+1;
                    break;
                }
            }
            gmsh::line ln;
            if( ln_tag == -1 )
            {
                ln_tag = new_ln_tag;
                ln.tag() = new_ln_tag++;
                ln.start_point_tag() = start_pt_tag;
                ln.end_point_tag() = end_pt_tag;
                lines.push_back(ln);
                map_firstp2ln[start_pt_tag].push_back(ln.tag()-1);
            }
            if( reverse ) sf.line_tags().push_back(-ln_tag);
            else          sf.line_tags().push_back(ln_tag);
        }
        surfaces.push_back(sf);
    }
    std::vector<gmsh::volume> volumes;
    for( int i=0; i<meshes.size(); i++ )
    {
        gmsh::volume vol;
        vol.tag() = i+1;
        for( const auto &sf_tag : face_indices_in_volumes[i] )
        {
            vol.surface_tags().push_back(sf_tag+1);
        }
        volumes.push_back(vol);
    }
    //
    geo.set_points( points );
    geo.set_lines( lines );
    geo.set_surfaces( surfaces );
    geo.set_volumes( volumes );

    //+++
    /*
    std::cout << "points: " << std::endl;
    for( const auto &pt : points )
    {
        pt.output();
    }
    std::cout << "lines: " << std::endl;
    for( const auto &ln : lines )
    {
        ln.output();
    }
    std::cout << "surfaces: " << std::endl;
    for( const auto &sf : surfaces )
    {
        sf.output();
    }
    std::cout << "volumes: " << std::endl;
    for( const auto &vol : volumes )
    {
        vol.output();
    }
    */
    //---

    return geo;
}