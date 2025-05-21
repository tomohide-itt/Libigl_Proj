#include <igl/cotmatrix.h>
#include <Eigen/Dense>
#include <iostream>
#include <filesystem>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_face_normals.h>
#include <igl/jet.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/boost/graph/helpers.h>
#include "gmsh.h"
#include "mesh.h"

namespace PMP = CGAL::Polygon_mesh_processing;

using CGAL_Point = CGAL::Simple_cartesian<double>::Point_3;
using CGAL_Mesh = CGAL::Surface_mesh<CGAL_Point>;

void integrate_mesh( const std::vector<migl::mesh>& meshes, const std::vector<Eigen::MatrixXd>& colors,
                     migl::mesh& all_mesh, Eigen::MatrixXd& all_colors );

void convert_to_cgal_mesh( const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, CGAL_Mesh& mesh );

void convert_to_eigen_mesh( const CGAL_Mesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F );

void detect_constraints( const CGAL_Mesh& mesh, std::map<CGAL_Mesh::Edge_index,bool>& constraints, double sharp_edge_angle_deg = 30.0 );

void append_mesh( const CGAL_Mesh& src_mesh, CGAL_Mesh& merged_mesh );

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

  //　すべての地盤メッシュに対して，頂点の最大値，最小値を表示する
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
  std::vector<migl::mesh> ground_meshes_out_remesh(ground_meshes.size());
  std::vector<Eigen::MatrixXd> colors_out_remesh(ground_meshes.size());
  {
    //int i = 0;
    for( int i=0; i<ground_meshes.size(); i++ )
    {
      std::cout << "remeshing ground_meshes_out[" << i << "] / " << ground_meshes_out.size()-1 << std::endl;
      // Ligigl -> CGAL
      CGAL_Mesh mesh;
      convert_to_cgal_mesh( ground_meshes_out[i].vertex_matrix(), ground_meshes_out[i].face_matrix(), mesh );
      std::cout << "  converted to CGAL mesh" << std::endl;
      // リメッシュ対象面の収集
      std::vector<CGAL_Mesh::Face_index> faces_to_remesh;
      for( auto f : mesh.faces() )
      {
        faces_to_remesh.push_back(f);
      }
      //リメッシュ実行
      double target_edge_length = 3;
      unsigned int iterations = 5;
      std::map<CGAL_Mesh::Edge_index,bool> constraints;
      detect_constraints( mesh, constraints, 1 );
      //PMP::isotropic_remeshing( faces_to_remesh, target_edge_length, mesh, PMP::parameters::number_of_iterations(iterations).protect_constraints(true) );
      PMP::isotropic_remeshing( faces_to_remesh, target_edge_length, mesh, PMP::parameters::number_of_iterations(iterations).edge_is_constrained_map(boost::make_assoc_property_map(constraints)) );
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
  //{
  //  std::cout << "remeshing all_ground_meshes_out" << std::endl;
  //  // Ligigl -> CGAL
  //  CGAL_Mesh mesh;
  //  convert_to_cgal_mesh( all_ground_mesh_out.vertex_matrix(), all_ground_mesh_out.face_matrix(), mesh );
  //  std::cout << "  converted to CGAL mesh" << std::endl;
  //  // リメッシュ対象面の収集
  //  std::vector<CGAL_Mesh::Face_index> faces_to_remesh;
  //  for( auto f : mesh.faces() )
  //  {
  //    faces_to_remesh.push_back(f);
  //  }
  //  //リメッシュ実行
  //  double target_edge_length = 3;
  //  unsigned int iterations = 5;
  //  std::map<CGAL_Mesh::Edge_index,bool> constraints;
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
    CGAL_Mesh merged_mesh;
    CGAL_Mesh mesh1, mesh13, mesh17;
    convert_to_cgal_mesh( ground_meshes_out[1].vertex_matrix(), ground_meshes_out[1].face_matrix(), mesh1 );
    convert_to_cgal_mesh( ground_meshes_out[13].vertex_matrix(), ground_meshes_out[13].face_matrix(), mesh13 );
    //convert_to_cgal_mesh( ground_meshes_out[17].vertex_matrix(), ground_meshes_out[17].face_matrix(), mesh17 );
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

void convert_to_cgal_mesh( const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, CGAL_Mesh& mesh )
{
  std::vector<CGAL_Mesh::Vertex_index> vertex_indices;
  for( int i=0; i<V.rows(); ++i )
  {
    vertex_indices.push_back( mesh.add_vertex( CGAL_Point(V(i,0), V(i,1), V(i,2))));
  }
  for( int i=0; i<F.rows(); ++i )
  {
    mesh.add_face( vertex_indices[F(i,0)], vertex_indices[F(i,1)], vertex_indices[F(i,2)] );
  }
}

void convert_to_eigen_mesh( const CGAL_Mesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F )
{
  std::map<CGAL_Mesh::Vertex_index, int> vertex_map;
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
    CGAL_Mesh::Halfedge_index h = mesh.halfedge(f);
    F(fi,0) = vertex_map[mesh.target(h)];
    h = mesh.next(h);
    F(fi,1) = vertex_map[mesh.target(h)];
    h = mesh.next(h);
    F(fi,2) = vertex_map[mesh.target(h)];
    ++fi;
  }
}

void detect_constraints( const CGAL_Mesh& mesh, std::map<CGAL_Mesh::Edge_index,bool>& constraints, double sharp_edge_angle_deg )
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
    if( f1 == CGAL_Mesh::null_face() || f2 == CGAL_Mesh::null_face()) continue;

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

/*
void append_mesh( const CGAL_Mesh& src_mesh, CGAL_Mesh& merged_mesh )
{
  std::map<CGAL_Mesh::Vertex_index, CGAL_Mesh::Vertex_index> vmap;
  for( auto v : src_mesh.vertices() )
  {
    auto p = src_mesh.point(v);
    //　既存の同一位置の頂点があれば使う（許容誤差付き）
    CGAL_Mesh::Vertex_index mv;
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
    std::vector<CGAL_Mesh::Vertex_index> verts;
    for( auto v : CGAL::vertices_around_face( src_mesh.halfedge(f), src_mesh ) )
    {
      verts.push_back(vmap[v]);
    }
    // 頂点数チェック
    std::set<CGAL_Mesh::Vertex_index> unique_verts(verts.begin(), verts.end());
    if (unique_verts.size() < 3) {
      std::cerr << "⚠️ 面スキップ：重複頂点による無効面" << std::endl;
      continue;
    }

    // 面追加
    auto face = merged_mesh.add_face(verts);
    if (face == CGAL_Mesh::null_face()) {
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