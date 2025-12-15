#include <igl/cotmatrix.h>
#include <Eigen/Dense>
#include <iostream>
#include <filesystem>
#include <cfloat>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_face_normals.h>
#include <igl/jet.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include "gmsh.h"
#include "mesh.h"
#include "cgal.h"

//namespace PMP = CGAL::Polygon_mesh_processing;

//using Point = CGAL::Simple_cartesian<double>::Point_3;
//using Vector = CGAL::Simple_cartesian<double>::Vector_3;
//using Mesh = CGAL::Surface_mesh<Point>;

void integrate_mesh( const std::vector<migl::mesh>& meshes, const std::vector<Eigen::MatrixXd>& colors,
                     migl::mesh& all_mesh, Eigen::MatrixXd& all_colors );

void convert_to_Mesh( const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Mesh& mesh );

void convert_to_eigen_mesh( const Mesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F );

void detect_constraint_edges( const Mesh& mesh, std::map<Mesh::Edge_index,bool>& constraints, double sharp_edge_angle_deg = 30.0 );
void detect_border_corners( const Mesh& mesh, std::map<Mesh::Vertex_index,bool>& constraints, double sharp_angle_deg = 30.0 );

void append_mesh( const Mesh& src_mesh, Mesh& merged_mesh );

bool extract_contact_surface( const Mesh &mesh1, const Mesh &mesh2, std::vector<Mesh> &contact_surfaces );
bool extract_contact_surface2( const Mesh &mesh1, const Mesh &mesh2, std::vector<Mesh> &contact_surfaces );
bool extract_contact_surfaces( const migl::mesh &mesh1, const migl::mesh &mesh2, std::vector<migl::mesh> &contact_surfaces );

bool extract_contact_edge( const Mesh &surf1, const Mesh &surf2, std::vector<std::vector<Mesh::Edge_index>> &contact_edges1, std::vector<std::vector<Mesh::Edge_index>> &contact_edges2 );
bool extract_contact_point_pair( const Mesh &surf1, const Mesh &surf2,
  std::vector<std::vector<std::pair<Mesh::Point,Mesh::Point>>> &contact_ppair1,
  std::vector<std::vector<std::pair<Mesh::Point,Mesh::Point>>> &contact_ppair2 );


void insert_contact_edge( const std::vector<std::pair<Mesh::Point,Mesh::Point>> &contact_ppair, const Mesh &rmesh, Mesh &mesh );

template<typename Mesh>
void remove_isolated_vertices(Mesh& mesh) {
    std::set<typename Mesh::Vertex_index> used_vertices;
    for (auto f : mesh.faces()) {
        for (auto v : CGAL::vertices_around_face(mesh.halfedge(f), mesh)) {
            used_vertices.insert(v);
        }
    }

    for (auto v : mesh.vertices()) {
        if (used_vertices.find(v) == used_vertices.end()) {
            mesh.remove_vertex(v);
        }
    }

    mesh.collect_garbage(); // 必須！
}

//=============================================================================================================
// メイン関数
int main(int argc, char *argv[])
{
  // アウトプットファイルをおくディレクトリを作成
  std::filesystem::path output_dir = std::filesystem::current_path() / "output";
  if (!std::filesystem::exists(output_dir))
  {
    std::filesystem::create_directory(output_dir);
  }
  //デバッグログ用のディレクトリを作成
  std::filesystem::path debug_dir = DEBUG_PATH;
  if (!std::filesystem::exists(debug_dir))
  {
    std::filesystem::create_directory(debug_dir);
  }
  //
  gmsh::geo ground_geo; //地盤モデル用geoオブジェクト
  gmsh::geo box_geo; //ボックス用geoオブジェクト
  // geoファイルを読み込む
  ground_geo.read(DATA_PATH "/ground_model.geo");
  box_geo.read(DATA_PATH "/box_shift_x6m.geo");
  // メッシュに変換する
  std::vector<migl::mesh> ground_meshes = ground_geo.get_meshes();
  std::vector<migl::mesh> box_meshes = box_geo.get_meshes();
  // すべてのメッシュに対して，面を外向きに整列させる 
  for( auto &mesh : ground_meshes) mesh.orient_faces_outward();
  for( auto &mesh : box_meshes)    mesh.orient_faces_outward();

  //　すべての地盤メッシュで，頂点の最小値を計算する
  Eigen::Vector3d min_vertex;
  {
    min_vertex << DBL_MAX, DBL_MAX, DBL_MAX;
    for( int i=0; i<ground_meshes.size(); i++ )
    {
      Eigen::Vector3d minv = ground_meshes[i].min_vertex();
      for( int k=0; k<3; k++ )
      {
        if( min_vertex(k) > minv(k) ) min_vertex(k) = minv(k);
      }
    }
    std::cout << "min point : " << min_vertex.transpose() << std::endl;
  }

  //地盤メッシュとボックスメッシュの頂点座標をmin_vertex分だけシフトする
  {
    for( auto &mesh : ground_meshes ) mesh.shift_vertices( -min_vertex );
    for( auto &mesh : box_meshes ) mesh.shift_vertices( -min_vertex );
  }

  //+++
  /*
  for( int i = 0; i < ground_meshes.size(); ++i )
  {
    std::cout << "ground_meshes[" << i << "] max: " << ground_meshes[i].max_vertex().transpose() << std::endl;
    std::cout << "ground_meshes[" << i << "] min: " << ground_meshes[i].min_vertex().transpose() << std::endl;
  }
  */
  //---

  //　すべての地盤メッシュに対して，V, Fを表示する
  //+++
  /*
  for( int i = 0; i < ground_meshes.size(); ++i )
  {
    std::cout << "ground_meshes[" << i << "] V: " << std::endl << ground_meshes[i].vertex_matrix() << std::endl;
    std::cout << "ground_meshes[" << i << "] F: " << std::endl << ground_meshes[i].face_matrix() << std::endl;
  }
  */
  //----
  //　地盤のメッシュを統合してgeoオブジェクトに変換し，geoファイルに出力
  //gmsh::geo integrated_ground_geo = migl::mesh::get_gmsh_geo(ground_meshes);
  //std::string integrated_ground_geo_filename = (output_dir / "integrated_ground_model.geo").string();
  //integrated_ground_geo.write(integrated_ground_geo_filename);

  //　すべての地盤メッシュの色を設定する
  std::vector<Eigen::MatrixXd> colors(ground_meshes.size());
  for( int i = 0; i < ground_meshes.size(); ++i )
  {
    Eigen::VectorXd I = Eigen::VectorXd::Constant(ground_meshes[i].vertex_matrix().rows(), static_cast<double>(i));
    igl::jet(I,0,ground_meshes.size(),colors[i]);
  }

  //　すべての地盤メッシュを統一する
  migl::mesh all_ground_mesh;
  Eigen::MatrixXd all_colors;
  integrate_mesh(ground_meshes, colors, all_ground_mesh, all_colors);

  //　すべての地盤メッシュとボックスメッシュとのブーリアン演算を行う
  std::vector<migl::mesh> ground_meshes_out(ground_meshes.size());
  std::vector<Eigen::MatrixXd> colors_out(ground_meshes.size());
  for( int i = 0; i < ground_meshes.size(); ++i )
  {
    igl::copyleft::cgal::mesh_boolean(ground_meshes[i].vertex_matrix(), ground_meshes[i].face_matrix(),
                            box_meshes[0].vertex_matrix(), box_meshes[0].face_matrix(),
                            igl::MESH_BOOLEAN_TYPE_MINUS,
                            ground_meshes_out[i].vertex_matrix(), ground_meshes_out[i].face_matrix());
    Eigen::VectorXd I = Eigen::VectorXd::Constant(ground_meshes_out[i].vertex_matrix().rows(), static_cast<double>(i));
    igl::jet(I,0,ground_meshes.size(),colors_out[i]);
  }

  // 各地盤メッシュをそれぞれgeoファイルに出力する 
  {
    for( int i=0; i<ground_meshes_out.size(); i++ )
    {
      std::vector<migl::mesh> partial_ground_meshes(1);
      partial_ground_meshes[0] = ground_meshes_out[i];
      gmsh::geo partial_ground_geo = migl::mesh::get_gmsh_geo(partial_ground_meshes);
      std::ostringstream ss_fname;
      std::string strfname = (output_dir / "ground_out_").string();
      ss_fname << strfname << i << ".geo";
      partial_ground_geo.write(ss_fname.str());
    }
  }

  //　地盤のメッシュを統合してgeoオブジェクトに変換し，geoファイルに出力する
  //gmsh::geo integrated_ground_out_geo = migl::mesh::get_gmsh_geo(ground_meshes_out);
  //std::string integrated_ground_out_geo_filename = (output_dir / "integrated_ground_out_model.geo").string();
  //integrated_ground_out_geo.write(integrated_ground_out_geo_filename);

  /*
  //　ブーリアン演算後の地盤メッシュの一部を統合してgeoファイルに出力する
  {
    std::vector<migl::mesh> partial_ground_meshes(1);
    partial_ground_meshes[0] = ground_meshes_out[5];
    //partial_ground_meshes[1] = ground_meshes_out[10];
    //partial_ground_meshes[2] = ground_meshes_out[12];
    gmsh::geo integrated_p_ground_geo = migl::mesh::get_gmsh_geo(partial_ground_meshes);
    std::string integrated_p_ground_geo_filename = (output_dir / "integrated_partial_ground_out_model.geo").string();
    integrated_p_ground_geo.write(integrated_p_ground_geo_filename);
  }
  */

  //CGALのisotropic_remeshingを使って，ブーリアン演算後の地盤メッシュをリメッシュ
  //時間がかかるのでとりあえず，コメントアウト
  /*
  std::vector<migl::mesh> ground_meshes_out_remesh(ground_meshes.size());
  std::vector<Eigen::MatrixXd> colors_out_remesh(ground_meshes.size());
  {
    //int i = 0;
    for( int i=0; i<ground_meshes.size(); i++ )
    {
      std::cout << "remeshing ground_meshes_out[" << i << "] / " << ground_meshes_out.size()-1 << std::endl;
      // Ligigl -> CGAL
      Mesh mesh;
      convert_to_Mesh( ground_meshes_out[i].vertex_matrix(), ground_meshes_out[i].face_matrix(), mesh );
      std::cout << "  converted to CGAL mesh" << std::endl;
      // リメッシュ対象面の収集
      std::vector<Mesh::Face_index> faces_to_remesh;
      for( auto f : mesh.faces() )
      {
        faces_to_remesh.push_back(f);
      }
      //リメッシュ実行
      double target_edge_length = 3;
      unsigned int iterations = 5;
      std::map<Mesh::Edge_index,bool> constraints;
      detect_constraint_edges( mesh, constraints, 1 );
      //PMP::isotropic_remeshing( faces_to_remesh, target_edge_length, mesh, PMP::parameters::number_of_iterations(iterations).protect_constraints(true) );
      PMP::isotropic_remeshing( faces_to_remesh, target_edge_length, mesh,
        PMP::parameters::number_of_iterations(iterations)
        .edge_is_constrained_map(boost::make_assoc_property_map(constraints))
        //.protect_constraints(true)
      );
      std::cout << "  remeshed" << std::endl;
      // CGAL -> Libigl
      convert_to_eigen_mesh( mesh, ground_meshes_out_remesh[i].vertex_matrix(), ground_meshes_out_remesh[i].face_matrix() );
      std::cout << "  converted to Libigl mesh" << std::endl;
      std::cout << "  finished" << std::endl;
    }
    //色の決定
    for( int j=0; j<ground_meshes.size(); j++ )
    {
      Eigen::VectorXd I = Eigen::VectorXd::Constant(ground_meshes_out_remesh[j].vertex_matrix().rows(), static_cast<double>(j));
      igl::jet(I,0,ground_meshes.size(),colors_out_remesh[j]);
    }
  }

  //　ボックスとのブーリアン演算を行ったすべての地盤メッシュを統一する
  migl::mesh all_ground_mesh_out;
  Eigen::MatrixXd all_colors_out;
  integrate_mesh(ground_meshes_out, colors_out, all_ground_mesh_out, all_colors_out);

  migl::mesh all_ground_mesh_out_remesh;
  Eigen::MatrixXd all_colors_out_remesh;
  integrate_mesh(ground_meshes_out_remesh, colors_out_remesh, all_ground_mesh_out_remesh, all_colors_out_remesh );
  */

  //{
  //  std::cout << "remeshing all_ground_meshes_out" << std::endl;
  //  // Ligigl -> CGAL
  //  Mesh mesh;
  //  convert_to_Mesh( all_ground_mesh_out.vertex_matrix(), all_ground_mesh_out.face_matrix(), mesh );
  //  std::cout << "  converted to CGAL mesh" << std::endl;
  //  // リメッシュ対象面の収集
  //  std::vector<Mesh::Face_index> faces_to_remesh;
  //  for( auto f : mesh.faces() )
  //  {
  //    faces_to_remesh.push_back(f);
  //  }
  //  //リメッシュ実行
  //  double target_edge_length = 3;
  //  unsigned int iterations = 5;
  //  std::map<Mesh::Edge_index,bool> constraints;
  //  detect_constraints( mesh, constraints, 10 );
  //  PMP::isotropic_remeshing( faces_to_remesh, target_edge_length, mesh, PMP::parameters::number_of_iterations(iterations).edge_is_constrained_map(boost::make_assoc_property_map(constraints)) );
  //  std::cout << "  remeshed" << std::endl;
  //  // CGAL -> Libigl
  //  convert_to_eigen_mesh( mesh, all_ground_mesh_out_remesh.vertex_matrix(), all_ground_mesh_out_remesh.face_matrix() );
  //  std::cout << "  converted to Libigl mesh" << std::endl;
  //  std::cout << "  finished" << std::endl;
  //  //色の決定
  //  Eigen::VectorXd I = Eigen::VectorXd::Constant(all_ground_mesh_out_remesh.vertex_matrix().rows(), static_cast<double>(1));
  //  igl::jet(I,0,ground_meshes.size(),all_colors_out_remesh);
  //}

  //メッシュの統合をGCALで試す →　append_meshで，接合部分の面が追加できない（既存のジオメトリと不整合であるため）
  /*
  {
    Mesh merged_mesh;
    Mesh mesh1, mesh13, mesh17;
    convert_to_Mesh( ground_meshes_out[1].vertex_matrix(), ground_meshes_out[1].face_matrix(), mesh1 );
    convert_to_Mesh( ground_meshes_out[13].vertex_matrix(), ground_meshes_out[13].face_matrix(), mesh13 );
    //convert_to_Mesh( ground_meshes_out[17].vertex_matrix(), ground_meshes_out[17].face_matrix(), mesh17 );
    append_mesh(mesh1,merged_mesh);
    append_mesh(mesh13,merged_mesh);
    //append_mesh(mesh17,merged_mesh);
    std::vector<migl::mesh> igl_merged_mesh(1);
    convert_to_eigen_mesh( merged_mesh, igl_merged_mesh[0].vertex_matrix(), igl_merged_mesh[0].face_matrix() );
    igl::opengl::glfw::Viewer viewer2;
    viewer2.data().set_mesh(igl_merged_mesh[0].vertex_matrix(), igl_merged_mesh[0].face_matrix() );
    viewer2.launch();
  }
  */
  //リメッシュをmesh1とmesh13とmesh17で試行
  {
    /*
    int num_mesh = 3;
    std::vector<migl::mesh> ground_meshes( num_mesh );
    ground_meshes[0] = ground_meshes_out[ 1];
    ground_meshes[1] = ground_meshes_out[13];
    ground_meshes[2] = ground_meshes_out[17];

    //各メッシュ間の接面を抽出 
    std::vector<migl::mesh> contact_surfaces;
    for( int i=0; i<num_mesh; i++ )
    {
      for( int j=i+1; j<num_mesh; j++ )
      {
        std::cout << "extracting contact surfaces between ground_meshes[" << i << "] and [" << j << "]" << std::endl;  
        extract_contact_surfaces( ground_meshes[i], ground_meshes[j], contact_surfaces );
        std::cout << "  finished" << std::endl;
      }
    }

    //contact_surfaceをリメッシュ target_edge_length= OK: 1, 1.5, NG: 2
    std::vector<migl::mesh> remeshed_contact_surfaces( contact_surfaces.size() );
    for( int i=0; i<contact_surfaces.size(); i++ )
    {
      std::cout << "remeshing contact_surfaces[" << i << "]" << std::endl;

      Mesh mesh;
      convert_to_Mesh( contact_surfaces[i].V(), contact_surfaces[i].F(), mesh );

      double target_edge_length = 1.5;
      unsigned int iterations = 2;
      std::map<Mesh::Edge_index,bool> constraint_edges;
      detect_constraint_edges( mesh, constraint_edges, 1 );
      std::map<Mesh::Vertex_index,bool> constraint_vertices;
      detect_border_corners( mesh, constraint_vertices, 1 );

      //+++
      if( i == 0 )
      {
        for( const auto &[e,b] : constraint_edges )
        {
          auto h = mesh.halfedge(e);
          auto v1 = mesh.source(h);
          auto v2 = mesh.target(h);
          auto p1 = mesh.point(v1);
          auto p2 = mesh.point(v2);
          std::cout << "constraints edge :" << e << " " << v1 << " ( " << p1 << " ) - " << v2 << " ( " << p2 << " )" << std::endl;  
        }
        for( const auto &[v,b] : constraint_vertices )
        {
          auto p = mesh.point(v);
          std::cout << "constraints vertex : " << v << " ( " << p << " )" << std::endl;
        }
      }
      //---

      PMP::isotropic_remeshing( faces(mesh), target_edge_length, mesh,
        PMP::parameters::number_of_iterations(iterations)
        .edge_is_constrained_map(boost::make_assoc_property_map(constraint_edges))
        .vertex_is_constrained_map(boost::make_assoc_property_map(constraint_vertices))
      );

      convert_to_eigen_mesh( mesh, remeshed_contact_surfaces[i].V(), remeshed_contact_surfaces[i].F() );

      std::cout << "  finished" << std::endl;
    }

    //ビューアー表示 
    {
      igl::opengl::glfw::Viewer viewer2;
      int id=0;
      bool remeshed=false;
      const auto &update_viewer = [&]()
      {
        viewer2.data().clear();
        if( !remeshed )
        {
          viewer2.data().set_mesh(contact_surfaces[id].V(), contact_surfaces[id].F() );
        }
        else
        {
          viewer2.data().set_mesh(remeshed_contact_surfaces[id].vertex_matrix(), remeshed_contact_surfaces[id].face_matrix() );
        }
      };
      update_viewer();
      viewer2.callback_key_pressed = [&](igl::opengl::glfw::Viewer &viewer2, unsigned char key, int modifiers)
      {
        if (key == '[')      id = (id+1) % contact_surfaces.size();
        else if (key == ']') id = (id-1+contact_surfaces.size()) % contact_surfaces.size();
        else if( key == 'r' && remeshed ) remeshed = false;
        else if( key == 'r' && !remeshed ) remeshed = true;
        else return false;
        update_viewer();
        return true;
      };
      viewer2.launch();
    }
    */

    int num_mesh = 3;
    std::vector<migl::mesh> ground_meshes( num_mesh );
    //for( int i=0; i<num_mesh; i++ ) ground_meshes[i] = ground_meshes_out[i];
    ground_meshes[0] = ground_meshes_out[ 1];
    ground_meshes[1] = ground_meshes_out[13];
    ground_meshes[2] = ground_meshes_out[17];
    //
    //CGALメッシュへ変換 
    std::vector<Mesh> cmeshes(num_mesh);
    for( int i=0; i<num_mesh; i++ ) convert_to_Mesh( ground_meshes[i].V(), ground_meshes[i].F(), cmeshes[i] );

    //各メッシュ間の接面を抽出して contact_surfaces に格納する
    std::vector<Mesh> contact_surfaces;
    std::vector<std::pair<int,int>> contact_mesh_indices;
    for( int i=0; i<num_mesh; i++ )
    {
      for( int j=i+1; j<num_mesh; j++ )
      {
        std::cout << "extracting contact surfaces between ground_meshes[" << i << "] and [" << j << "]" << std::endl;
        Mesh surface = cgal::extract_contact_surface( cmeshes[i], cmeshes[j] );
        if( surface.num_faces() != 0 )
        {
          contact_surfaces.push_back( surface );
          contact_mesh_indices.push_back( std::make_pair(i,j) );
        }
        std::cout << "  finished" << std::endl;
      }
    }
    std::cout << "done extracting contact surfaces" << std::endl;

    // contact_surfaces を Eigen meshに変換
    std::vector<migl::mesh> igl_contact_surfaces( contact_surfaces.size() );
    for( int i=0; i<contact_surfaces.size(); i++ )
    {
        convert_to_eigen_mesh( contact_surfaces[i], igl_contact_surfaces[i].V(), igl_contact_surfaces[i].F() );
    }

    //+++ 各 contact_sufaces　が　どのメッシュとどのメッシュの接合面かを出力 
    std::cout << "contact_mesh_indices" << std::endl;
    for( int i=0; i<contact_surfaces.size(); i++ )
    {
      std::cout << "contact_mesh_indices[" << i << "] = " << "(" << contact_mesh_indices[i].first << ", " << contact_mesh_indices[i].second << ")" << std::endl;
    }
    //---

    //contact_surfaces間で共通する境界線を持つものを特定し， neigh_contact_surface_indices に格納する 
    std::vector<std::vector<int>> neigh_contact_surface_indices( contact_surfaces.size() );
    //std::vector<std::vector<std::vector<Mesh::Edge_index>>> contact_edges( contact_surfaces.size() );
    std::vector<std::vector<std::vector<std::pair<Mesh::Point,Mesh::Point>>>> contact_point_pair( contact_surfaces.size() );
    for( int i=0; i<contact_surfaces.size(); i++ )
    {
      for( int j=i+1; j<contact_surfaces.size(); j++ )
      {
        std::cout << "extracting contact edges between contact_surfaces[" << i << "] and [" << j << "]" << std::endl; 
        //if( extract_contact_edge( contact_surfaces[i], contact_surfaces[j], contact_edges[i], contact_edges[j] ) )
        //{
        //  neigh_contact_surface_indices[i].push_back(j);
        //  neigh_contact_surface_indices[j].push_back(i);
        //  std::cout << "  found contact edges" << std::endl;
        //}
        //else
        //{
        //  std::cout << "  not found contact edges" << std::endl;
        //}
        if( extract_contact_point_pair( contact_surfaces[i], contact_surfaces[j], contact_point_pair[i], contact_point_pair[j]) )
        {
          neigh_contact_surface_indices[i].push_back(j);
          neigh_contact_surface_indices[j].push_back(i);
          std::cout << "  found contact ppair" << std::endl;
        }
        else
        {
          std::cout << "  not found contact ppair" << std::endl;
        }
      }
    }
    std::cout << "done extracting contact ppair" << std::endl;
    //+++
    std::cout << "neigh_contact_surface_indices" << std::endl;
    for( int i=0; i<contact_surfaces.size(); i++ )
    {
      std::cout << "neigh_contact_surface_indices[" << i << "] = ( ";
      for( const auto &neigh_contact_surface_index : neigh_contact_surface_indices[i] )
      {
        std::cout << neigh_contact_surface_index << " ";
      }
      std::cout << ")" << std::endl;
    }
    std::cout << "number of contact_edges" << std::endl;
    for( int i=0; i<contact_surfaces.size(); i++ )
    {
      std::cout << "number of contact_edges[" << i << "] = ( ";
      for( const auto &ppairs : contact_point_pair[i] )
      {
        std::cout << ppairs.size() << " ";
      }
      std::cout << ")" << std::endl;
    }
    //---
    //contact_surfaceをリメッシュ target_edge_length= OK: 1, 1.5, NG: 2
    double target_edge_len = 1.5;
    std::cout << "target_edge_length >> ";
    std::cin >> target_edge_len;
    std::vector<bool> remeshed( contact_surfaces.size(), false );
    for( int i=0; i<contact_surfaces.size(); i++ )
    {
      std::cout << "remeshing contact_surfaces[" << i << "]" << std::endl;

      //リメッシュ実行
      double target_edge_length = target_edge_len;
      unsigned int iterations = 2;
      std::map<Mesh::Edge_index,bool> constraint_edges;
      detect_constraint_edges( contact_surfaces[i], constraint_edges, 1 );
      std::map<Mesh::Vertex_index,bool> constraint_vertices;
      detect_border_corners( contact_surfaces[i], constraint_vertices, 1 );

      PMP::isotropic_remeshing( faces(contact_surfaces[i]), target_edge_length, contact_surfaces[i],
        PMP::parameters::number_of_iterations(iterations)
        .edge_is_constrained_map(boost::make_assoc_property_map(constraint_edges))
        .vertex_is_constrained_map(boost::make_assoc_property_map(constraint_vertices))
      );
      std::cout << "  done remeshing" << std::endl;
      remeshed[i] = true;
      //
      //nei_contact_surface_indices[i]の中に既にリメッシュされたcontact_surfaceがあるか調べ，
      //あれば，そのリメッシュされたcontact_surfaceの contact_edge上の点をコピー
      for( int j=0; j<neigh_contact_surface_indices[i].size(); j++ )
      {
        int index = neigh_contact_surface_indices[i][j];
        if( remeshed[index] )
        {
          std::cout << "inserting contact edge contact_surfaces[" << index << "] -> [" << i << "]" << std::endl;
          insert_contact_edge( contact_point_pair[i][j], contact_surfaces[index], contact_surfaces[i] );
        }
      }
    }
    //リメッシュされたcontact surfacesを Eigen meshに変換
    std::vector<migl::mesh> igl_remeshed_contact_surfaces( contact_surfaces.size() );
    for( int i=0; i<contact_surfaces.size(); i++ )
    {
        convert_to_eigen_mesh( contact_surfaces[i], igl_remeshed_contact_surfaces[i].V(), igl_remeshed_contact_surfaces[i].F() );
    }
    //---
    //ビューアー表示 
    {
      igl::opengl::glfw::Viewer viewer2;
      int id=0;
      bool remeshed=false;
      const auto &update_viewer = [&]()
      {
        viewer2.data().clear();
        if( !remeshed )
        {
          viewer2.data().set_mesh(igl_contact_surfaces[id].vertex_matrix(), igl_contact_surfaces[id].face_matrix() );
        }
        else
        {
          viewer2.data().set_mesh(igl_remeshed_contact_surfaces[id].vertex_matrix(), igl_remeshed_contact_surfaces[id].face_matrix() );
        }
      };
      update_viewer();
      viewer2.callback_key_pressed = [&](igl::opengl::glfw::Viewer &viewer2, unsigned char key, int modifiers)
      {
        if (key == '[')      id = (id+1) % igl_contact_surfaces.size();
        else if (key == ']') id = (id-1+igl_contact_surfaces.size()) % igl_contact_surfaces.size();
        else if( key == 'r' && remeshed ) remeshed = false;
        else if( key == 'r' && !remeshed ) remeshed = true;
        else return false;
        update_viewer();
        return true;
      };
      viewer2.launch();
    }
    //テスト □ はOK，□ □ はNG　→　接合面は離れた独立した面となることもあるため，対応が必要 
    {
      migl::mesh migl_mesh(3,4,2);
      migl_mesh.vertex_matrix() << 0, 0, 0,
                                   1, 0, 0,
                                   1, 1, 0,
                                   0, 1, 0;
                                   //2, 0, 0,
                                   //3, 0, 0,
                                   //3, 1, 0,
                                   //2, 1, 0;
      migl_mesh.face_matrix() << 0, 1, 2,
                                 0, 2, 3;
                                 //4, 5, 6,
                                 //4, 6, 7;
      Mesh Mesh;
      convert_to_Mesh( migl_mesh.vertex_matrix(), migl_mesh.face_matrix(), Mesh );
      PMP::isotropic_remeshing( faces(Mesh), 0.1, Mesh,
        PMP::parameters::number_of_iterations(2)
        //.edge_is_constrained_map(boost::make_assoc_property_map(constraints))
      );
      migl::mesh migl_remesh;
      convert_to_eigen_mesh( Mesh, migl_remesh.vertex_matrix(), migl_remesh.face_matrix() );
      igl::opengl::glfw::Viewer viewer3;
      viewer3.data().set_mesh(migl_remesh.vertex_matrix(), migl_remesh.face_matrix() );
      viewer3.launch();
    }
  }

  /*

  // ビューアの初期化
  igl::opengl::glfw::Viewer viewer;
  int mesh_id = 0;
  bool all_meshes = false;
  bool boolean_done = false;
  bool remesh = false;
  //
  // ビューアの情報を更新するラムダ関数
  const auto &update = [&]()
  {
    viewer.data().clear();
    if( all_meshes )
    {
      if( boolean_done )
      {
        if( remesh )
        {
          viewer.data().set_mesh(all_ground_mesh_out_remesh.vertex_matrix(), all_ground_mesh_out_remesh.face_matrix());
          viewer.data().set_colors(all_colors_out_remesh);
        }
        else
        {
          viewer.data().set_mesh(all_ground_mesh_out.vertex_matrix(), all_ground_mesh_out.face_matrix());
          viewer.data().set_colors(all_colors_out);
        }
      }
      else
      {
        viewer.data().set_mesh(all_ground_mesh.vertex_matrix(), all_ground_mesh.face_matrix());
        viewer.data().set_colors(all_colors);
      }
    }
    else
    {
      if( boolean_done )
      {
        if( remesh )
        {
          viewer.data().set_mesh(ground_meshes_out_remesh[mesh_id].vertex_matrix(), ground_meshes_out_remesh[mesh_id].face_matrix());
          viewer.data().set_colors(colors_out_remesh[mesh_id]);
        }
        else
        {
          viewer.data().set_mesh(ground_meshes_out[mesh_id].vertex_matrix(), ground_meshes_out[mesh_id].face_matrix());
          viewer.data().set_colors(colors_out[mesh_id]);
        }
      }
      else
      {
        viewer.data().set_mesh(ground_meshes[mesh_id].vertex_matrix(), ground_meshes[mesh_id].face_matrix());
        viewer.data().set_colors(colors[mesh_id]);
      }
    }
  };
  update();

  // コールバック関数を登録する
  // (]キーでメッシュを前進, [キーでメッシュを後退，@キーですべてのメッシュを表示/非表示)
  viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifiers)
  {
    if (key == '[')      mesh_id = (mesh_id+1) % ground_meshes.size();
    else if (key == ']') mesh_id = (mesh_id-1+ground_meshes.size()) % ground_meshes.size();
    else if( key == '@' && all_meshes )  all_meshes = false;
    else if( key == '@' && !all_meshes ) all_meshes = true;
    else if( key == 'b' && boolean_done ) boolean_done = false;
    else if( key == 'b' && !boolean_done ) boolean_done = true;
    else if( key == 'r' && remesh ) remesh = false;
    else if( key == 'r' && !remesh ) remesh = true;
    else                 return false;
    update();
    return true;
  };
  viewer.launch();

  */

  return 0;
}

void integrate_mesh( const std::vector<migl::mesh>& meshes, const std::vector<Eigen::MatrixXd>& colors,
                     migl::mesh& all_mesh, Eigen::MatrixXd& all_colors )
{
  int V_rows = 0;
  int F_rows = 0;
  for( int i = 0; i < meshes.size(); ++i )
  {
    V_rows += meshes[i].vertex_matrix().rows();
    F_rows += meshes[i].face_matrix().rows();
  }
  all_mesh.vertex_matrix().resize(V_rows, 3);
  all_mesh.face_matrix().resize(F_rows, 3);
  all_colors.resize(V_rows, 3);
  int v_offset = 0;
  int f_offset = 0;
  for( int i = 0; i < meshes.size(); ++i )
  {
    for( int j = 0; j < meshes[i].vertex_matrix().rows(); ++j )
    {
      all_mesh.vertex_matrix().row(v_offset+j) = meshes[i].vertex_matrix().row(j);
      all_colors.row(v_offset+j) = colors[i].row(j);
    }
    for( int j = 0; j < meshes[i].face_matrix().rows(); ++j )
    {
      for( int k=0; k<3; k++ )
      {
        all_mesh.face_matrix()(f_offset+j,k) = meshes[i].face_matrix()(j,k) + v_offset;
      }
    }
    v_offset += meshes[i].vertex_matrix().rows();
    f_offset += meshes[i].face_matrix().rows();
  }
}

void convert_to_Mesh( const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Mesh& mesh )
{
  std::vector<Mesh::Vertex_index> vertex_indices;
  for( int i=0; i<V.rows(); ++i )
  {
    vertex_indices.push_back( mesh.add_vertex( Point(V(i,0), V(i,1), V(i,2))));
  }
  for( int i=0; i<F.rows(); ++i )
  {
    mesh.add_face( vertex_indices[F(i,0)], vertex_indices[F(i,1)], vertex_indices[F(i,2)] );
  }
}

void convert_to_eigen_mesh( const Mesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F )
{
  std::map<Mesh::Vertex_index, int> vertex_map;
  V.resize( mesh.number_of_vertices(), 3);
  int vi = 0;
  for( auto v : mesh.vertices() )
  {
    const auto& p = mesh.point(v);
    V.row(vi) = Eigen::RowVector3d( p.x(), p.y(), p.z() );
    vertex_map[v] = vi++;
  }
  F.resize( mesh.number_of_faces(), 3);
  int fi = 0;
  for( auto f : mesh.faces() )
  {
    Mesh::Halfedge_index h = mesh.halfedge(f);
    F(fi,0) = vertex_map[mesh.target(h)];
    h = mesh.next(h);
    F(fi,1) = vertex_map[mesh.target(h)];
    h = mesh.next(h);
    F(fi,2) = vertex_map[mesh.target(h)];
    ++fi;
  }
}

void detect_constraint_edges( const Mesh& mesh, std::map<Mesh::Edge_index,bool>& constraints, double sharp_edge_angle_deg )
{
  //
  for( auto e : mesh.edges() )
  {
    if( mesh.is_border(e) )
    {
      constraints[e] = true;
      continue;
    }
    //
    auto h = mesh.halfedge(e);
    auto f1 = mesh.face(h);
    auto f2 = mesh.face(mesh.opposite(h));
    if( f1 == Mesh::null_face() || f2 == Mesh::null_face()) continue;

    auto n1 = PMP::compute_face_normal(f1, mesh);
    auto n2 = PMP::compute_face_normal(f2, mesh);

    double dot = n1 * n2;
    double norm1 = std::sqrt(n1.squared_length());
    double norm2 = std::sqrt(n2.squared_length());
    double cos_angle = dot / (norm1 * norm2);

    if( cos_angle < std::cos(sharp_edge_angle_deg * CGAL_PI / 180.0 ) )
    {
      constraints[e] = true;
    }
  }
}

void detect_border_corners( const Mesh& mesh, std::map<Mesh::Vertex_index,bool>& constraints, double sharp_angle_deg )
{
  std::vector<Mesh::Edge_index> borders;
  for( const auto &e : mesh.edges() )
  {
    if( mesh.is_border(e) ) borders.push_back(e);
  }

  for( const auto &border : borders )
  {
    auto h1 = mesh.halfedge(border);
    auto v1s = mesh.source(h1);
    auto v1t = mesh.target(h1);
    auto p1s = mesh.point(v1s);
    auto p1t = mesh.point(v1t);
    auto vec1 = p1t - p1s;
    Vector vec2;
    for( const auto &border : borders )
    {
      auto h = mesh.halfedge(border);
      auto vs = mesh.source(h);
      auto ps = mesh.point(vs);
      if( CGAL::squared_distance(p1t, ps) < 1e-10 )
      {
        auto vt = mesh.target(h);
        auto pt = mesh.point(vt);
        vec2 = pt - ps;
        break;
      }
    }

    double dot = vec1 * vec2;
    double norm1 = std::sqrt(vec1.squared_length());
    double norm2 = std::sqrt(vec2.squared_length());
    double cos_angle = dot / (norm1 * norm2);

    if( cos_angle < std::cos(sharp_angle_deg * CGAL_PI / 180.0 ) )
    {
      constraints[v1t] = true;
    }

  }
}

/*
void append_mesh( const Mesh& src_mesh, Mesh& merged_mesh )
{
  std::map<Mesh::Vertex_index, Mesh::Vertex_index> vmap;
  for( auto v : src_mesh.vertices() )
  {
    auto p = src_mesh.point(v);
    //　既存の同一位置の頂点があれば使う（許容誤差付き）
    Mesh::Vertex_index mv;
    bool found = false;
    for( auto u : merged_mesh.vertices() )
    {
      if( CGAL::squared_distance(p, merged_mesh.point(u)) < 1e-10 )
      {
        mv = u;
        found = true;
        break;
      }
    }
    if( !found ) mv = merged_mesh.add_vertex(p);
    vmap[v] = mv;
  }
  for( auto f : src_mesh.faces() )
  {
    std::vector<Mesh::Vertex_index> verts;
    for( auto v : CGAL::vertices_around_face( src_mesh.halfedge(f), src_mesh ) )
    {
      verts.push_back(vmap[v]);
    }
    // 頂点数チェック
    std::set<Mesh::Vertex_index> unique_verts(verts.begin(), verts.end());
    if (unique_verts.size() < 3) {
      std::cerr << "⚠️ 面スキップ：重複頂点による無効面" << std::endl;
      continue;
    }

    // 面追加
    auto face = merged_mesh.add_face(verts);
    if (face == Mesh::null_face()) {
      std::cerr << "⚠️ 面追加失敗：";
      for (auto v : verts) std::cerr << v << " ";
      std::cerr << std::endl;
      for( auto v : verts )
      {
        std::cout << v << "(" << merged_mesh.point(v) << ") ";
      }
      std::cout << std::endl;
      Eigen::Vector3d p0 = Eigen::Vector3d(merged_mesh.point(verts[0]).x(), merged_mesh.point(verts[0]).y(), merged_mesh.point(verts[0]).z() );
      Eigen::Vector3d p1 = Eigen::Vector3d(merged_mesh.point(verts[1]).x(), merged_mesh.point(verts[1]).y(), merged_mesh.point(verts[1]).z() );
      Eigen::Vector3d p2 = Eigen::Vector3d(merged_mesh.point(verts[2]).x(), merged_mesh.point(verts[2]).y(), merged_mesh.point(verts[2]).z() );

      double area = 0.5 * ((p1 - p0).cross(p2 - p0)).norm();
      std::cout << "Area = " << area << std::endl;
    }
  }
}
*/

bool extract_contact_surface( const Mesh &mesh1, const Mesh &mesh2, std::vector<Mesh> &contact_surfaces )
{
  double eps = 1e-10;
  Mesh surface;
  std::map<Mesh::Vertex_index, Mesh::Vertex_index> vmap;
  for( auto f1 : mesh1.faces() )
  {
    std::vector<Point> pnts1;
    std::vector<Mesh::Vertex_index> verts1;
    for( auto v : CGAL::vertices_around_face( mesh1.halfedge(f1), mesh1 ) )
    {
      pnts1.push_back(mesh1.point(v));
      verts1.push_back(v);
    }
    //
    for( auto f2 : mesh2.faces() )
    {
      std::vector<Point> pnts2;
      std::vector<Mesh::Vertex_index> verts2;
      for( auto v : CGAL::vertices_around_face( mesh2.halfedge(f2), mesh2 ) )
      {
        pnts2.push_back(mesh2.point(v));
        verts2.push_back(v);
      }
      //+++
      //std::cout << "face1: " << f1 << " num = " << mesh1.num_faces() << std::endl;
      //std::cout << "  " << verts1[0] << ": " << mesh1.point(verts1[0]) << std::endl;
      //std::cout << "  " << verts1[1] << ": " << mesh1.point(verts1[1]) << std::endl;
      //std::cout << "  " << verts1[2] << ": " << mesh1.point(verts1[2]) << std::endl;
      //std::cout << "face2: " << f2 << " num = " << mesh2.num_faces() << std::endl;
      //std::cout << "  " << verts2[0] << ": " << mesh2.point(verts2[0]) << std::endl;
      //std::cout << "  " << verts2[1] << ": " << mesh2.point(verts2[1]) << std::endl;
      //std::cout << "  " << verts2[2] << ": " << mesh2.point(verts2[2]) << std::endl;
      //---
      if( ( CGAL::squared_distance(pnts1[0], pnts2[0]) < eps && CGAL::squared_distance(pnts1[1], pnts2[1]) < eps && CGAL::squared_distance(pnts1[2], pnts2[2]) < eps ) ||
          ( CGAL::squared_distance(pnts1[0], pnts2[0]) < eps && CGAL::squared_distance(pnts1[1], pnts2[2]) < eps && CGAL::squared_distance(pnts1[2], pnts2[1]) < eps ) ||
          ( CGAL::squared_distance(pnts1[0], pnts2[1]) < eps && CGAL::squared_distance(pnts1[1], pnts2[0]) < eps && CGAL::squared_distance(pnts1[2], pnts2[2]) < eps ) ||
          ( CGAL::squared_distance(pnts1[0], pnts2[1]) < eps && CGAL::squared_distance(pnts1[1], pnts2[2]) < eps && CGAL::squared_distance(pnts1[2], pnts2[0]) < eps ) ||
          ( CGAL::squared_distance(pnts1[0], pnts2[2]) < eps && CGAL::squared_distance(pnts1[1], pnts2[0]) < eps && CGAL::squared_distance(pnts1[2], pnts2[1]) < eps ) ||
          ( CGAL::squared_distance(pnts1[0], pnts2[2]) < eps && CGAL::squared_distance(pnts1[1], pnts2[1]) < eps && CGAL::squared_distance(pnts1[2], pnts2[0]) < eps ) )
      {
        //+++
        //std::cout << "same triangle " << std::endl;
        //---
        std::vector<Mesh::Vertex_index> verts(3);
        for( int i=0; i<3; i++ )
        {
          if( vmap.count(verts1[i]) ) verts[i] = vmap[verts1[i]];
          else
          {
            verts[i] = surface.add_vertex(pnts1[i]);
            vmap[verts1[i]] = verts[i];
          }
        }
        //+++
        //std::cout << "verts: " << verts[0] << ", " << verts[1] << ", " << verts[2] << std::endl;
        //std::cout << "  " << verts[0] << ": " << surface.point(verts[0]) << std::endl;
        //std::cout << "  " << verts[1] << ": " << surface.point(verts[1]) << std::endl;
        //std::cout << "  " << verts[2] << ": " << surface.point(verts[2]) << std::endl;
        //---
        surface.add_face(verts);
        break;
      }
    }
  }
  //+++
  //std::cout << "surface.num_faces()=" << surface.num_faces() << std::endl; 
  //---
  if( surface.num_faces() != 0 )
  {
    // 修復処理
    PMP::triangulate_faces(surface);
    PMP::remove_degenerate_faces(surface);
    remove_isolated_vertices(surface);
    PMP::duplicate_non_manifold_vertices(surface);

    if( !CGAL::is_valid_polygon_mesh(surface) )
    {
      std::cerr << "Generated surface is invalid!" << std::endl;
    }
    for (auto h : halfedges(surface)) {
      bool is_border_h = CGAL::is_border(h, surface);
      bool is_border_opposite = CGAL::is_border(opposite(h, surface), surface);
      if (is_border_h && is_border_opposite) {
          std::cerr << "Broken halfedge: both sides are border! h=" << h << std::endl;
      }
    }

    contact_surfaces.push_back(surface);
    return true;
  }
  return false;
}

bool extract_contact_surface2( const Mesh &mesh1, const Mesh &mesh2, std::vector<Mesh> &contact_surfaces )
{
  Mesh surface;
  Mesh mesh1_copy = mesh1;
  Mesh mesh2_copy = mesh2;
  try {
    PMP::corefine_and_compute_intersection(mesh1_copy, mesh2_copy, surface);
  } catch (const std::exception& e) {
    std::cerr << "Corefine failed: " << e.what() << std::endl;
    return false;
  }

  if (surface.is_empty() || surface.num_faces() == 0) return false;

  // optional: 後処理
  PMP::triangulate_faces(surface);
  PMP::remove_degenerate_faces(surface);
  PMP::remove_isolated_vertices(surface);
  PMP::duplicate_non_manifold_vertices(surface);

  if (!CGAL::is_valid_polygon_mesh(surface)) {
    std::cerr << "Extracted surface is invalid!" << std::endl;
    return false;
  }

  contact_surfaces.push_back(surface);
  return true;
}

bool extract_contact_surfaces( const migl::mesh &mesh1, const migl::mesh &mesh2, std::vector<migl::mesh> &contact_surfaces )
{
  //+++
  //std::cout << "extract_contact_surfacecs" << std::endl;
  //std::cout << "  mesh1.num_face: " << mesh1.F().rows();
  //std::cout << "  mesh2.num_face: " << mesh2.F().rows();
  //---
  double eps = 1e-10;
  std::map<int, int> vmap;
  std::vector<Eigen::RowVector3d> verts;
  std::vector<Eigen::Vector3i> faces; 
  //
  for( int fi1=0; fi1<mesh1.F().rows(); fi1++ )
  {
    std::vector<int> vindices1;
    std::vector<Eigen::RowVector3d> pnts1;
    for( int k=0; k<3; k++ ) vindices1.push_back( mesh1.F(fi1,k) );
    for( int k=0; k<3; k++ ) pnts1.push_back( mesh1.V(vindices1[k]) );
    //
    for( int fi2=0; fi2<mesh2.F().rows(); fi2++ )
    {
      std::vector<int> vindices2;
      std::vector<Eigen::RowVector3d> pnts2;
      for( int k=0; k<3; k++ ) vindices2.push_back( mesh2.F(fi2,k) );
      for( int k=0; k<3; k++ ) pnts2.push_back( mesh2.V(vindices2[k]) );

      if( ( pnts1[0].isApprox(pnts2[0],eps) && pnts1[1].isApprox(pnts2[1],eps) && pnts1[2].isApprox(pnts2[2],eps) ) ||
          ( pnts1[0].isApprox(pnts2[0],eps) && pnts1[1].isApprox(pnts2[2],eps) && pnts1[2].isApprox(pnts2[1],eps) ) ||
          ( pnts1[0].isApprox(pnts2[1],eps) && pnts1[1].isApprox(pnts2[0],eps) && pnts1[2].isApprox(pnts2[2],eps) ) ||
          ( pnts1[0].isApprox(pnts2[1],eps) && pnts1[1].isApprox(pnts2[2],eps) && pnts1[2].isApprox(pnts2[0],eps) ) ||
          ( pnts1[0].isApprox(pnts2[2],eps) && pnts1[1].isApprox(pnts2[0],eps) && pnts1[2].isApprox(pnts2[1],eps) ) ||
          ( pnts1[0].isApprox(pnts2[2],eps) && pnts1[1].isApprox(pnts2[1],eps) && pnts1[2].isApprox(pnts2[0],eps) ) )
      {
        //+++
        //std::cout << "  fi1 = " << fi1 << ", fi2 = " << fi2 << std::endl;
        //---
        Eigen::Vector3i face = Eigen::Vector3i::Zero(3);
        for( int k=0; k<3; k++ )
        {
          if( vmap.count( vindices1[k] ) )
          {
            face(k) = vmap[vindices1[k]];
          }
          else
          {
            int vindex = verts.size();
            verts.push_back( pnts1[k] );
            face(k) = vindex;
            vmap[vindices1[k]] = vindex;
          }
        }
        //+++
        //std::cout << "  face: " << face << std::endl;
        //---
        faces.push_back( face );
        break;
      }
    }
  }
  if( faces.size() != 0 )
  {
    migl::mesh surface( 3, verts.size(), faces.size() );
    for( int i=0; i<verts.size(); i++ )
    {
      surface.V(i,0) = verts[i](0);
      surface.V(i,1) = verts[i](1);
      surface.V(i,2) = verts[i](2);
    }
    for( int i=0; i<faces.size(); i++ )
    {
      surface.F(i,0) = faces[i](0);
      surface.F(i,1) = faces[i](1);
      surface.F(i,2) = faces[i](2);
    }
    contact_surfaces.push_back( surface );
    return true;
  }
  return false;
}

bool extract_contact_edge( const Mesh &surf1, const Mesh &surf2,
  std::vector<std::vector<Mesh::Edge_index>> &contact_edges1,
  std::vector<std::vector<Mesh::Edge_index>> &contact_edges2 )
{
  double eps = 1e-10;
  std::vector<Mesh::Edge_index> contact_edges;
  for( const auto &edge1 : surf1.edges() )
  {
    auto h1 = surf1.halfedge(edge1);
    Mesh::Vertex_index v11 = surf1.source(h1);
    Mesh::Vertex_index v12 = surf1.target(h1);
    Mesh::Point p11 = surf1.point(v11);
    Mesh::Point p12 = surf1.point(v12);
    for( const auto &edge2 : surf2.edges() )
    {
      auto h2 = surf2.halfedge(edge2);
      Mesh::Vertex_index v21 = surf2.source(h2);
      Mesh::Vertex_index v22 = surf2.target(h2);
      Mesh::Point p21 = surf2.point(v21);
      Mesh::Point p22 = surf2.point(v22);
      if( ( CGAL::squared_distance(p11, p21) < eps && CGAL::squared_distance(p12, p22) < eps ) ||
          ( CGAL::squared_distance(p11, p22) < eps && CGAL::squared_distance(p12, p21) < eps ) )
      {
        contact_edges.push_back(edge1);
        break;
      }
    }
  }
  if( contact_edges.size() != 0 )
  {
    contact_edges1.push_back( contact_edges );
    contact_edges2.push_back( contact_edges );
    return true;
  }
  return false;
}

bool extract_contact_point_pair( const Mesh &surf1, const Mesh &surf2,
  std::vector<std::vector<std::pair<Mesh::Point,Mesh::Point>>> &contact_ppair1,
  std::vector<std::vector<std::pair<Mesh::Point,Mesh::Point>>> &contact_ppair2 )
{
  double eps = 1e-10;
  std::vector<std::pair<Mesh::Point,Mesh::Point>> contact_ppair;
  for( const auto &edge1 : surf1.edges() )
  {
    auto h1 = surf1.halfedge(edge1);
    Mesh::Vertex_index v11 = surf1.source(h1);
    Mesh::Vertex_index v12 = surf1.target(h1);
    Mesh::Point p11 = surf1.point(v11);
    Mesh::Point p12 = surf1.point(v12);
    for( const auto &edge2 : surf2.edges() )
    {
      auto h2 = surf2.halfedge(edge2);
      Mesh::Vertex_index v21 = surf2.source(h2);
      Mesh::Vertex_index v22 = surf2.target(h2);
      Mesh::Point p21 = surf2.point(v21);
      Mesh::Point p22 = surf2.point(v22);
      if( ( CGAL::squared_distance(p11, p21) < eps && CGAL::squared_distance(p12, p22) < eps ) ||
          ( CGAL::squared_distance(p11, p22) < eps && CGAL::squared_distance(p12, p21) < eps ) )
      {
        contact_ppair.push_back(std::make_pair(p11,p12));
        break;
      }
    }
  }
  if( contact_ppair.size() != 0 )
  {
    contact_ppair1.push_back( contact_ppair );
    contact_ppair2.push_back( contact_ppair );
    return true;
  }
  return false;
}

void insert_contact_edge( const std::vector<std::pair<Mesh::Point,Mesh::Point>> &contact_ppair, const Mesh &src_mesh, Mesh &dist_mesh )
{
  double eps = 1e-10;
  //src_meshのborderのうちcontact_ppair（無向エッジの集合）に含まれる点を収集する
  std::vector<Mesh::Point> p_on_contact_edges_src;
  for( const auto &v : src_mesh.vertices() )
  {
    if( !src_mesh.is_border(v) ) continue;
    auto p = src_mesh.point(v);
    for( const auto &ppair : contact_ppair )
    {
      auto p1 = ppair.first;
      auto p2 = ppair.second;
      //両端のどちらかと一致
      if( CGAL::squared_distance(p1, p) < eps || CGAL::squared_distance(p2, p) < eps )
      {
        p_on_contact_edges_src.push_back(p);
        break;
      }
      //辺上にある 
      auto vec21 = p2 - p1;
      auto vecp1 = p - p1;
      double dot = vec21 * vecp1;
      double norm1 = std::sqrt(vec21.squared_length());
      double norm2 = std::sqrt(vecp1.squared_length());
      double cos = dot / (norm1 * norm2);
      if( fabs(cos - 1.0) < eps && norm1 > norm2 )
      {
        p_on_contact_edges_src.push_back(p);
        break;
      }
    }
  }
  //+++
  std::cout << "p_on_contact_edges_src: " << std::endl;
  for( const auto & p : p_on_contact_edges_src )
  {
    std::cout << "  ( " << p << " )" << std::endl;
  }
  //---
  //dist_meshのborderのうちcontact_ppair（無向エッジの集合）に含まれる点を収集する
  std::vector<Mesh::Point> p_on_contact_edges_dist;
  for( const auto &v : dist_mesh.vertices() )
  {
    if( !dist_mesh.is_border(v) ) continue;
    auto p = dist_mesh.point(v);
    for( const auto &ppair : contact_ppair )
    {
      auto p1 = ppair.first;
      auto p2 = ppair.second;
      //両端のどちらかと一致
      if( CGAL::squared_distance(p1, p) < eps || CGAL::squared_distance(p2, p) < eps )
      {
        p_on_contact_edges_dist.push_back(p);
        break;
      }
      //辺上にある 
      auto vec21 = p2 - p1;
      auto vecp1 = p - p1;
      double dot = vec21 * vecp1;
      double norm1 = std::sqrt(vec21.squared_length());
      double norm2 = std::sqrt(vecp1.squared_length());
      double cos = dot / (norm1 * norm2);
      if( fabs(cos - 1.0) < eps && norm1 > norm2 )
      {
        p_on_contact_edges_dist.push_back(p);
        break;
      }
    }
  }
  //+++
  std::cout << "p_on_contact_edges_dist: " << std::endl;
  for( const auto & p : p_on_contact_edges_dist )
  {
    std::cout << "  ( " << p << " )" << std::endl;
  }
  //---

  //p_on_contact_edges_srcのうち，p_on_contact_distに含まれないものを収集する
  std::vector<Mesh::Point> pnts_add;
  for( const auto &p_src : p_on_contact_edges_src )
  {
    bool found = false;
    for( const auto &p_dist : p_on_contact_edges_dist )
    {
      if( CGAL::squared_distance(p_src, p_dist) < eps )
      {
        found = true;
        break;
      }
    }
    if( !found )
    {
      pnts_add.push_back( p_src );
    }
  }
  //+++
  std::cout << "pnts_add: " << std::endl;
  for( const auto & p : pnts_add )
  {
    std::cout << "  ( " << p << " )" << std::endl;
  }
  //---

  //p_on_contact_edges_distのうち，p_on_contact_srcに含まれないものを収集する
  std::vector<Mesh::Point> pnts_delete;
  for( const auto &p_dist : p_on_contact_edges_dist )
  {
    bool found = false;
    for( const auto &p_src : p_on_contact_edges_src )
    {
      if( CGAL::squared_distance(p_src, p_dist) < eps )
      {
        found = true;
        break;
      }
    }
    if( !found )
    {
      pnts_delete.push_back( p_dist );
    }
  }
  //+++
  std::cout << "pnts_delete: " << std::endl;
  for( const auto & p : pnts_delete )
  {
    std::cout << "  ( " << p << " )" << std::endl;
  }
  //---

  if( pnts_add.size() == 0 && pnts_delete.size() == 0 ) return;

  //dist_meshのborderにあるhalfedgeを収集する
  std::vector<Mesh::Halfedge_index> border_halfedges;
  for( const auto &h : dist_mesh.halfedges() )
  {
    if( dist_mesh.is_border(h) ) border_halfedges.push_back(h);
  }

  //p_addをdist_meshに追加し，新しく面を作る 
  if( pnts_add.size() != 0 )
  {
    for( auto &h : border_halfedges )
    {
      Mesh::Vertex_index v1 = dist_mesh.source(h);
      Mesh::Vertex_index v2 = dist_mesh.target(h);
      auto p1 = dist_mesh.point(v1);
      auto p2 = dist_mesh.point(v2);
      for( const auto &p : pnts_add )
      {
        //辺上にあるかcheck 
        auto vec21 = p2 - p1;
        auto vecp1 = p - p1;
        double dot = vec21 * vecp1;
        double norm1 = std::sqrt(vec21.squared_length());
        double norm2 = std::sqrt(vecp1.squared_length());
        double cos = dot / (norm1 * norm2);
        //この辺上になければcontinue
        if( fabs(cos - 1.0) >= eps || norm1 <= norm2 ) continue;
        //hを含む面で，v1でもv2でもない頂点を特定
        Mesh::Vertex_index v3;
        for( const auto &v : vertices_around_face(h,dist_mesh) )
        {
          if( v != v1 && v != v2 )
          {
            v3 = v; break;
          }
        }
        //面を消去して新しい2つの面を作る
        Mesh::Vertex_index vp = dist_mesh.add_vertex(p);
        CGAL::Euler::remove_face( h, dist_mesh );
        CGAL::Euler::add_face( std::vector<Mesh::Vertex_index>{vp, v2, v3}, dist_mesh );
        CGAL::Euler::add_face( std::vector<Mesh::Vertex_index>{vp, v3, v1}, dist_mesh );
        //+++
        std::cout << "point : ( " << p << " ) is added to dist_mesh" << std::endl;
        //---
        break;
      }
    }
  }

  //pnts_deleteをdist_meshから削除して面を統合する 
  if( pnts_delete.size() != 0 )
  {

  }
}
