#include <igl/cotmatrix.h>
#include <Eigen/Dense>
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_face_normals.h>
#include <igl/jet.h>
#include "libs/gmsh/gmsh.h"
#include "mesh.h"

int main(int argc, char *argv[])
{
  gmsh::geo geo;
  // geoファイルを読み込む
  geo.read(DATA_PATH "/ground_model.geo");
  // Volume毎にメッシュを変換してベクトル化する
  std::vector<migl::mesh> meshes = geo.get_meshes();
  // すべてのメッシュに対して，面を外向きに整列させる 
  for( auto &mesh : meshes) mesh.orient_faces_outward();
  
  igl::opengl::glfw::Viewer viewer;
  int mesh_id = 0;
  //ビューアの情報を更新するラムダ関数
  const auto &update = [&]()
  {
    viewer.data().clear();
    viewer.data().set_mesh(meshes[mesh_id].vertex_matrix(), meshes[mesh_id].face_matrix());
    Eigen::VectorXd I = Eigen::VectorXd::Constant(meshes[mesh_id].vertex_matrix().rows(), static_cast<double>(mesh_id));
    Eigen::MatrixXd C;
    igl::jet(I,0,meshes.size(),C);
    viewer.data().set_colors(C);
  };
  update();

  // コールバック関数を登録する
  // (]キーでメッシュを前進, [キーでメッシュを後退)
  viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifiers)
  {
    if (key == '[')      mesh_id = (mesh_id+1) % meshes.size();
    else if (key == ']') mesh_id = (mesh_id-1+meshes.size()) % meshes.size();
    else                 return false;
    update();
    return true;
  };
  viewer.launch();

  return 0;
}
