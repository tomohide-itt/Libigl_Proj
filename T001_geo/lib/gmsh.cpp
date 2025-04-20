#include "gmsh.h"

//read the geo file
void gmsh::geo::read( const std::string& filename )
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    // Read the file line by line
    // and store each line in a vector
    // (or process it as needed)
    std::string line;
    std::vector<std::string> lines;
    while (std::getline(file, line))
    {
        lines.push_back(line);
    }
    for (const auto& l : lines)
    {
        std::cout << l << std::endl;
    }
}