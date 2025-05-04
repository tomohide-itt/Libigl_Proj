#include <igl/cotmatrix.h>
#include <Eigen/Dense>
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_face_normals.h>
#include <igl/jet.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include "libs/gmsh/gmsh.h"
#include "mesh.h"

void cal_mesh( const std::vector<migl::mesh>& gmeshes, const migl::mesh& bmesh, const int mesh_id, const bool all_meshes,
               Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& Colors );

//=============================================================================================================
// メイン関数
int main(int argc, char *argv[])
{
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
  
  igl::opengl::glfw::Viewer viewer;
  int mesh_id = 0;
  bool all_meshes = false;
  //
  //ビューアの情報を更新するラムダ関数
  const auto &update = [&]()
  {
    viewer.data().clear();
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd Colors;
    cal_mesh( ground_meshes, box_meshes[0], mesh_id, all_meshes, V, F, Colors );
    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(Colors);
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
    else                 return false;
    update();
    return true;
  };
  viewer.launch();

  return 0;
}

void cal_mesh( const std::vector<migl::mesh>& gmeshes, const migl::mesh& bmesh, const int mesh_id, const bool all_meshes,
               Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& Colors )
{
  if( !all_meshes)
  {
    igl::copyleft::cgal::mesh_boolean(gmeshes[mesh_id].vertex_matrix(), gmeshes[mesh_id].face_matrix(),
                            bmesh.vertex_matrix(), bmesh.face_matrix(),
                            igl::MESH_BOOLEAN_TYPE_MINUS,
                            V, F);
    Eigen::VectorXd I = Eigen::VectorXd::Constant(V.rows(), static_cast<double>(mesh_id));
    igl::jet(I,0,gmeshes.size(),Colors);
  }
  else
  {
    std::vector<Eigen::MatrixXd> V_temp(gmeshes.size());
    std::vector<Eigen::MatrixXi> F_temp(gmeshes.size());
    std::vector<Eigen::MatrixXd> Colors_temp(gmeshes.size());
    int V_rows = 0;
    int F_rows = 0;
    for( int i = 0; i < gmeshes.size(); ++i )
    {
      igl::copyleft::cgal::mesh_boolean(gmeshes[i].vertex_matrix(), gmeshes[i].face_matrix(),
                              bmesh.vertex_matrix(), bmesh.face_matrix(),
                              igl::MESH_BOOLEAN_TYPE_MINUS,
                              V_temp[i], F_temp[i]);
      Eigen::VectorXd I = Eigen::VectorXd::Constant(V_temp[i].rows(), static_cast<double>(i));
      igl::jet(I,0,gmeshes.size(),Colors_temp[i]);
      V_rows += V_temp[i].rows();
      F_rows += F_temp[i].rows();
    }
    V.resize(V_rows, 3);
    F.resize(F_rows, 3);
    Colors.resize(V_rows, 3);
    int v_offset = 0;
    int f_offset = 0;
    for( int i = 0; i < gmeshes.size(); ++i )
    {
      for( int j = 0; j < V_temp[i].rows(); ++j )
      {
        V.row(v_offset+j) = V_temp[i].row(j);
        Colors.row(v_offset+j) = Colors_temp[i].row(j);
      }
      for( int j = 0; j < F_temp[i].rows(); ++j )
      {
        for( int k=0; k<3; k++ )
        {
          F(f_offset+j,k) = F_temp[i](j,k) + v_offset;
        }
      }
      v_offset += V_temp[i].rows();
      f_offset += F_temp[i].rows();
    }
  }
}
