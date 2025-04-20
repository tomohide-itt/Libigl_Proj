//#include <igl/readOFF.h>
//#include <igl/writeOBJ.h>
#include <iostream>
#include "lib/gmsh.h"

int main(int argc, char *argv[])
{

  std::cout << DATA_PATH "/ground_model.geo" << std::endl;
  gmsh::geo::read(DATA_PATH "/ground_model.geo");

  return 0;
}
