#include <igl/cotmatrix.h>
#include <Eigen/Dense>
#include <iostream>
#include <filesystem>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_face_normals.h>
#include <igl/jet.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include "gmsh.h"
#include "mesh.h"

void integrate_mesh( const std::vector<migl::mesh>& meshes, const std::vector<Eigen::MatrixXd>& colors,
                     migl::mesh& all_mesh, Eigen::MatrixXd& all_colors );

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
  box_geo.read(DATA_PATH "/box.geo");
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

  //　地盤のメッシュを統合してgeoオブジェクトに変換する
  gmsh::geo integrated_ground_geo = migl::mesh::get_gmsh_geo(ground_meshes);
  //　geoファイルを出力する
  std::string integrated_ground_geo_filename = (output_dir / "integrated_ground_model.geo").string();
  integrated_ground_geo.write(integrated_ground_geo_filename);

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

  //　地盤のメッシュを統合してgeoオブジェクトに変換する
  gmsh::geo integrated_ground_out_geo = migl::mesh::get_gmsh_geo(ground_meshes_out);
  //　geoファイルを出力する
  std::string integrated_ground_out_geo_filename = (output_dir / "integrated_ground_out_model.geo").string();
  integrated_ground_out_geo.write(integrated_ground_out_geo_filename);

  //　ボックスとのブーリアン演算を行ったすべての地盤メッシュを統一する
  migl::mesh all_ground_mesh_out;
  Eigen::MatrixXd all_colors_out;
  integrate_mesh(ground_meshes_out, colors_out, all_ground_mesh_out, all_colors_out);

  // ビューアの初期化
  igl::opengl::glfw::Viewer viewer;
  int mesh_id = 0;
  bool all_meshes = false;
  bool boolean_done = false;
  //
  // ビューアの情報を更新するラムダ関数
  const auto &update = [&]()
  {
    viewer.data().clear();
    if( all_meshes )
    {
      if( boolean_done)
      {
        viewer.data().set_mesh(all_ground_mesh_out.vertex_matrix(), all_ground_mesh_out.face_matrix());
        viewer.data().set_colors(all_colors_out);
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
        viewer.data().set_mesh(ground_meshes_out[mesh_id].vertex_matrix(), ground_meshes_out[mesh_id].face_matrix());
        viewer.data().set_colors(colors_out[mesh_id]);
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
