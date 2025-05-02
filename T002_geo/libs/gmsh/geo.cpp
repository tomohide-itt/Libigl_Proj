#include "gmsh.h"

// read the geo file
void gmsh::geo::read(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    std::cout << "Reading file: " << filename << std::endl;
    //
    // Calculate the number of points, lines, surfaces, and volumes
    // and the max point, line, surface, and volume tags
    read_max_tags(file);
    // Read the entities
    read_entities(file);
    // Close the file
    file.close();
    std::cout << "Finished reading file: " << filename << std::endl;
}


std::vector<migl::mesh> gmsh::geo::get_meshes() const
{
    std::vector<migl::mesh> meshes;

    // Create a mesh for each volume
    int dim = 3;
    for( const auto &volume : this->m_volumes )
    {
        int num_faces = volume.size();
        int num_vertices = 0;
        std::map<int, int> map_ptag2vindex;
        std::map<int, int> map_vindex2ptag;
        for( int i = 0; i < num_faces; i++ )
        {
            int sf_tag = volume.surface_tag(i);
            int sf_index = this->m_surface_map_tag2index.at(sf_tag);
            surface sf = this->m_surfaces[sf_index];
            int num_lines = sf.size();
            //+++
            //std::cout << "i: " << i << " sf_tag: " << sf_tag << " sf_index: " << sf_index << " num_lines: " << num_lines << std::endl;
            //---
            for( int j=0; j<num_lines; j++ )
            {
                int ln_tag = sf.line_tag(j);
                int ln_positive_tag = sf.line_positive_tag(j);
                //+++
                //std::cout << "  j: " << j << " ln_tag: " << ln_tag << std::endl;
                //---
                int ln_index = this->m_line_map_tag2index.at(ln_positive_tag);
                line ln = this->m_lines[ln_index];
                int pt1_index = this->m_point_map_tag2index.at(ln.start_point_tag());
                int pt2_index = this->m_point_map_tag2index.at(ln.end_point_tag());
                point pt1 = this->m_points[pt1_index];
                point pt2 = this->m_points[pt2_index];
                int pt1_tag = pt1.tag();
                int pt2_tag = pt2.tag();
                if( !map_ptag2vindex.count(pt1_tag) )
                {
                    map_ptag2vindex[pt1_tag] = num_vertices;
                    map_vindex2ptag[num_vertices] = pt1_tag;
                    num_vertices++;
                }
                if( !map_ptag2vindex.count(pt2_tag) )
                {
                    map_ptag2vindex[pt2_tag] = num_vertices;
                    map_vindex2ptag[num_vertices] = pt2_tag;
                    num_vertices++;
                }
            }
        }
        //+++
        //std::cout << "num_vertices: " << num_vertices << std::endl;
        //std::cout << "map_ptag2vindex: " << std::endl;
        //for( const auto &pair : map_ptag2vindex )
        //{
        //    std::cout << "  " << pair.first << " -> " << pair.second << std::endl;
        //}
        //std::cout << "map_vindex2ptag: " << std::endl;
        //for( const auto &pair : map_vindex2ptag )
        //{
        //    std::cout << "  " << pair.first << " -> " << pair.second << std::endl;
        //}
        //---
        // Create the vertex matrix
        migl::mesh mesh(dim, num_vertices, num_faces);
        for( int i=0; i<num_vertices; i++ )
        {
            int pt_tag = map_vindex2ptag.at(i);
            int pt_index = this->m_point_map_tag2index.at(pt_tag);
            point pt = this->m_points[pt_index];
            mesh.vertex_matrix()(i,0) = pt.x();
            mesh.vertex_matrix()(i,1) = pt.y();
            mesh.vertex_matrix()(i,2) = pt.z();
        }
        //+++
        std::cout << "mesh.vertex_matrix(): " << std::endl;
        std::cout << mesh.vertex_matrix() << std::endl;
        //---
        // Create the face matrix
        for( int i=0; i<num_faces; i++ )
        {
            int sf_tag = volume.surface_tag(i);
            int sf_index = this->m_surface_map_tag2index.at(sf_tag);
            surface sf = this->m_surfaces[sf_index];
            int num_lines = sf.size();
            //+++
            //std::cout << "i: " << i << " sf_tag: " << sf_tag << " sf_index: " << sf_index << " num_lines: " << num_lines << std::endl;
            //---
            for( int j=0; j<num_lines; j++ )
            {
                int ln_tag = sf.line_tag(j);
                int ln_positive_tag = sf.line_positive_tag(j);
                //+++
                //std::cout << "  j: " << j << " ln_tag: " << ln_tag << std::endl;
                //---
                int ln_index = this->m_line_map_tag2index.at(ln_positive_tag);
                line ln = this->m_lines[ln_index];
                int pt1_index = this->m_point_map_tag2index.at(ln.start_point_tag());
                int pt2_index = this->m_point_map_tag2index.at(ln.end_point_tag());
                point pt1 = this->m_points[pt1_index];
                point pt2 = this->m_points[pt2_index];
                int pt1_tag = pt1.tag();
                int pt2_tag = pt2.tag();
                //+++
                //std::cout << "    pt1_tag: " << pt1_tag << " pt2_tag: " << pt2_tag << std::endl;
                //---
                if( ln_tag > 0 )
                {
                    mesh.face_matrix()(i,j) = map_ptag2vindex.at(pt1_tag);
                    //+++
                    //std::cout << "    registered pt_tag: " << pt1_tag << std::endl;
                    //std::cout << "    registered vindex:" << map_ptag2vindex.at(pt1_tag) << std::endl;
                    //---
                }
                else
                {
                    mesh.face_matrix()(i,j) = map_ptag2vindex.at(pt2_tag);
                    //+++
                    //std::cout << "    registered pt_tag: " << pt2_tag << std::endl;
                    //std::cout << "    registered vindex:" << map_ptag2vindex.at(pt2_tag) << std::endl;
                    //---
                }
            }
        }
        //+++
        std::cout << "mesh.face_matrix(): " << std::endl;
        std::cout << mesh.face_matrix() << std::endl;
        //---
        meshes.push_back(mesh);
    }

    return meshes;
}

// read the tag of a point, line, surface, or volume
int gmsh::geo::read_tag(const std::string &str, const std::string &tag_name)
{
    std::string search_str = tag_name + "(";
    if (str.find(search_str) != std::string::npos)
    {
        size_t pos = str.find(search_str);
        size_t start_pos = str.find("(", pos);
        size_t end_pos = str.find(")", start_pos);
        std::string tag_str = str.substr(start_pos + 1, end_pos - start_pos - 1);
        int tag = std::stoi(tag_str);
        return tag;
    }
    return -1;
}

// read the max tags
void gmsh::geo::read_max_tags(std::ifstream &file)
{
    // Move the file pointer to the beginning
    file.seekg(0, std::ios::beg);
    // Initialize the number of points, lines, surfaces, and volumes
    this->m_num_points = 0;
    this->m_num_lines = 0;
    this->m_num_surfaces = 0;
    this->m_num_volumes = 0;
    this->m_max_point_tag = 0;
    this->m_max_line_tag = 0;
    this->m_max_surface_tag = 0;
    this->m_max_volume_tag = 0;
    //
    // Read the file line by line
    std::string str;
    while (std::getline(file, str))
    {
        // Check for the max point tag & the number of points
        if (str.find("Point(") != std::string::npos)
        {
            int tag = read_tag(str, "Point");
            if (tag > this->m_max_point_tag)
            {
                this->m_max_point_tag = tag;
            }
            this->m_num_points++;
        }
        // Check for the max line tag & the number of lines
        else if (str.find("Line(") != std::string::npos)
        {
            int tag = read_tag(str, "Line");
            if (tag > this->m_max_line_tag)
            {
                this->m_max_line_tag = tag;
            }
            this->m_num_lines++;
        }
        // Check for the max surface tag & the number of surfaces
        else if (str.find("Line Loop(") != std::string::npos)
        {
            int tag = read_tag(str, "Line Loop");
            if (tag > this->m_max_surface_tag)
            {
                this->m_max_surface_tag = tag;
            }
            this->m_num_surfaces++;
        }
        // Check for the max volume tag & the number of volumes
        else if (str.find("Surface Loop(") != std::string::npos)
        {
            int tag = read_tag(str, "Surface Loop");
            if (tag > this->m_max_volume_tag)
            {
                this->m_max_volume_tag = tag;
            }
            this->m_num_volumes++;
        }
    }
    // Output the max tags and the number of points
    std::cout << "  Max point tag: " << this->m_max_point_tag << std::endl;
    std::cout << "  Number of points: " << this->m_num_points << std::endl;
    std::cout << "  Max line tag: " << this->m_max_line_tag << std::endl;
    std::cout << "  Number of lines: " << this->m_num_lines << std::endl;
    std::cout << "  Max surface tag: " << this->m_max_surface_tag << std::endl;
    std::cout << "  Number of surfaces: " << this->m_num_surfaces << std::endl;
    std::cout << "  Max volume tag: " << this->m_max_volume_tag << std::endl;
    std::cout << "  Number of volumes: " << this->m_num_volumes << std::endl;
    // Move the file pointer to the beginning
    file.clear();
    file.seekg(0, std::ios::beg);
}

// read a point
void gmsh::geo::read_point(const std::string &str)
{
    // Get the point tag
    int tag = read_tag(str, "Point");
    // Get the point coordinates
    size_t pos = str.find("{");
    size_t end_pos = str.find(",", pos);
    std::string coords_str = str.substr(pos + 1, end_pos - pos - 1);
    double x = std::stod(coords_str);
    pos = end_pos + 1;
    end_pos = str.find(",", pos);
    std::string y_str = str.substr(pos, end_pos - pos);
    double y = std::stod(y_str);
    pos = end_pos + 1;
    end_pos = str.find("}", pos);
    std::string z_str = str.substr(pos, end_pos - pos);
    double z = std::stod(z_str);
    // Create a new point object and add it to the points vector
    int index = this->m_points.size();
    this->m_points.push_back(point(x, y, z, tag));
    // Add the point to the point map
    this->m_point_map_index2tag[index] = tag;
    this->m_point_map_tag2index[tag] = index;
}

// read a line
void gmsh::geo::read_line(const std::string &str)
{
    // Get the line tag
    int tag = read_tag(str, "Line");
    // Get the start and end point tags
    size_t pos = str.find("{");
    size_t end_pos = str.find(",", pos);
    std::string start_point_str = str.substr(pos + 1, end_pos - pos - 1);
    int start_point_tag = std::stoi(start_point_str);
    pos = end_pos + 1;
    end_pos = str.find("}", pos);
    std::string end_point_str = str.substr(pos, end_pos - pos);
    int end_point_tag = std::stoi(end_point_str);
    // Create a new line object and add it to the lines vector
    int index = this->m_lines.size();
    this->m_lines.push_back(line(tag, start_point_tag, end_point_tag));
    // Add the line to the line map
    this->m_line_map_index2tag[index] = tag;
    this->m_line_map_tag2index[tag] = index;
}

// read a surface
void gmsh::geo::read_surface(const std::string &str)
{
    // Get the surface tag
    int tag = read_tag(str, "Line Loop");
    // Get the line loop tags
    size_t pos = str.find("{");
    size_t end_pos = str.find(",", pos);
    std::string line_loop_str = str.substr(pos + 1, end_pos - pos - 1);
    std::vector<int> line_tags;
    while (end_pos != std::string::npos)
    {
        int line_tag = std::stoi(line_loop_str);
        line_tags.push_back(line_tag);
        pos = end_pos + 1;
        end_pos = str.find(",", pos);
        if (end_pos == std::string::npos)
        {
            end_pos = str.find("}", pos);
        }
        line_loop_str = str.substr(pos, end_pos - pos);
    }
    // Create a new surface object and add it to the surfaces vector
    int index = this->m_surfaces.size();
    this->m_surfaces.push_back(surface(tag, line_tags));
    // Add the surface to the surface map
    this->m_surface_map_index2tag[index] = tag;
    this->m_surface_map_tag2index[tag] = index;
}

// read a volume
void gmsh::geo::read_volume(const std::string &str)
{
    // Get the volume tag
    int tag = read_tag(str, "Surface Loop");
    // Get the surface loop tags
    size_t pos = str.find("{");
    size_t end_pos = str.find(",", pos);
    std::string surface_loop_str = str.substr(pos + 1, end_pos - pos - 1);
    std::vector<int> surface_tags;
    while (end_pos != std::string::npos)
    {
        int surface_tag = std::stoi(surface_loop_str);
        surface_tags.push_back(surface_tag);
        pos = end_pos + 1;
        end_pos = str.find(",", pos);
        if (end_pos == std::string::npos)
        {
            end_pos = str.find("}", pos);
        }
        surface_loop_str = str.substr(pos, end_pos - pos);
    }
    // Create a new volume object and add it to the volumes vector
    int index = this->m_volumes.size();
    this->m_volumes.push_back(volume(tag, surface_tags));
    // Add the volume to the volume map
    this->m_volume_map_index2tag[index] = tag;
    this->m_volume_map_tag2index[tag] = index;
}

// read the entities
void gmsh::geo::read_entities(std::ifstream &file)
{
    // Move the file pointer to the beginning
    file.clear();
    file.seekg(0, std::ios::beg);
    // Read the file line by line
    std::string str;
    while (std::getline(file, str))
    {
        // Check for the point tag
        if (str.find("Point(") != std::string::npos)
        {
            read_point(str);
        }
        // Check for the line tag
        else if (str.find("Line(") != std::string::npos)
        {
            read_line(str);
        }
        // Check for the surface tag
        else if (str.find("Line Loop(") != std::string::npos)
        {
            read_surface(str);
        }
        // Check for the volume tag
        else if (str.find("Surface Loop(") != std::string::npos)
        {
            read_volume(str);
        }
    }
    //+++
    /*
    {
        // Output the points, lines, surfaces, and volumes
        std::cout << "Points:" << std::endl;
        for (const auto &point : this->m_points)
        {
            point.output();
        }
        std::cout << "Lines:" << std::endl;
        for (const auto &line : this->m_lines)
        {
            line.output();
        }
        std::cout << "Surfaces:" << std::endl;
        for (const auto &surface : this->m_surfaces)
        {
            surface.output();
        }
        std::cout << "Volumes:" << std::endl;
        for (const auto &volume : this->m_volumes)
        {
            volume.output();
        }
        std::cout << "Finished reading entities." << std::endl;
    }
    */
    //---
}