#include <igl/cotmatrix.h>
#include <Eigen/Dense>
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_face_normals.h>
#include "libs/gmsh/gmsh.h"
#include "mesh.h"

int main(int argc, char *argv[])
{
  gmsh::geo geo;
  // Read the geo file
  geo.read(DATA_PATH "/ground_model.geo");
  //geo.read(DATA_PATH "/box2.geo");
  // Get the meshes from the geo object
  std::vector<migl::mesh> meshes = geo.get_meshes();


  meshes[0].orient_faces_outward();

  
  // Output the meshes
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(meshes[0].vertex_matrix(), meshes[0].face_matrix());
  viewer.data().set_face_based(true);
  viewer.launch();
  
  

  return 0;
}
