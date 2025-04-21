#include <igl/cotmatrix.h>
#include <Eigen/Dense>
#include <iostream>
#include "gmsh/gmsh.h"

int main(int argc, char *argv[])
{
  gmsh::geo geo;
  // Read the geo file
  geo.read(DATA_PATH "/ground_model.geo");

  

  return 0;
}
