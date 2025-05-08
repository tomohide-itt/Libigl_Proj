#ifndef GMESH_H
#define GMESH_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>
//#include "../migl/mesh.h"
//#include "mesh.h"

namespace migl
{
    class mesh; // Forward declaration of the mesh class
}

namespace gmsh
{
    class point
    {
        public:
            point() : m_x(0), m_y(0), m_z(0), m_tag(-1) {}
            point( double x, double y, double z, int tag );
            point( const point& p );
            point& operator=( const point& p );
            // Getters
            double x() const { return m_x; }
            double y() const { return m_y; }
            double z() const { return m_z; }
            int tag() const { return m_tag; }
            // Setters
            double& x() { return m_x; }
            double& y() { return m_y; }
            double& z() { return m_z; }
            int& tag() { return m_tag; }
            // Print function
            void output() const;

        private:
            double m_x;
            double m_y;
            double m_z;
            int m_tag;
    }; //point

    class line
    {
        public:
            line() : m_tag(-1) {}
            line( int tag, int start_point_tag, int end_point_tag );
            line( const line& l );
            line& operator=( const line& l );
            // Getters
            int tag() const { return m_tag; }
            int start_point_tag() const { return m_start_point_tag; }
            int end_point_tag() const { return m_end_point_tag; }
            // Setters
            int& tag() { return m_tag; }
            int& start_point_tag() { return m_start_point_tag; }
            int& end_point_tag() { return m_end_point_tag; }
            // Print function
            void output() const;

        private:
            int m_tag;
            int m_start_point_tag;
            int m_end_point_tag;
    }; //line

    class surface
    {
        public:
            surface() : m_tag(-1) {}
            surface( int tag, const std::vector<int>& line_tags );
            surface( const surface& s );
            surface& operator=( const surface& s );
            // Getters
            int tag() const { return m_tag; }
            int positive_tag() const;
            int size() const { return m_line_tags.size(); }
            int line_tag( int i ) const { return m_line_tags[i]; }
            int line_positive_tag( int i ) const { return m_line_tags[i] > 0 ? m_line_tags[i] : -m_line_tags[i]; }
            // Get the line tag without sign
            // Setters
            int& tag() { return m_tag; }
            int& line_tag( int i ) { return m_line_tags[i]; }
            std::vector<int>& line_tags() { return m_line_tags; }
            // Print function
            void output() const;

        private:
            int m_tag;
            std::vector<int> m_line_tags;
    }; //surface

    class volume
    {
        public:
            volume() : m_tag(-1) {}
            volume( int tag, const std::vector<int>& surface_tags );
            volume( const volume& v );
            volume& operator=( const volume& v );
            // Getters
            int tag() const { return m_tag; }
            int size() const { return m_surface_tags.size(); }
            int surface_tag( int i ) const { return m_surface_tags[i]; }
            // Setters
            int& tag() { return m_tag; }
            int& surface_tag( int i ) { return m_surface_tags[i]; }
            std::vector<int>& surface_tags() { return m_surface_tags; }
            // Print function
            void output() const;

        private:
            int m_tag;
            std::vector<int> m_surface_tags;
    }; //volume

    class geo
    {
        public:
            void read( const std::string& filename );
            void write( const std::string& filename );
            // Getters
            int num_points() const { return m_num_points; }
            int num_lines() const { return m_num_lines; }
            int num_surfaces() const { return m_num_surfaces; }
            int num_volumes() const { return m_num_volumes; }
            // Setters
            void set_points( const std::vector<point>& points ) { m_points = points; }
            void set_lines( const std::vector<line>& lines ) { m_lines = lines; }
            void set_surfaces( const std::vector<surface>& surfaces ) { m_surfaces = surfaces; }
            void set_volumes( const std::vector<volume>& volumes ) { m_volumes = volumes; }
            // Get the vector of meshes
            std::vector<migl::mesh> get_meshes() const;

        private:
            int read_tag( const std::string& str, const std::string& tag_name );
            void read_max_tags( std::ifstream& file );
            void read_point( const std::string& str );
            void read_line( const std::string& str );
            void read_surface( const std::string& str );
            void read_volume( const std::string& str );
            void read_entities( std::ifstream& file );
            void write_point( std::ofstream& file );
            void write_line( std::ofstream& file );
            void write_surface( std::ofstream& file );
            void write_volume( std::ofstream& file );
            void write_entities( std::ofstream& file );
        
        private:
            int m_max_point_tag = 0;
            int m_max_line_tag = 0;
            int m_max_surface_tag = 0;
            int m_max_volume_tag = 0;
            int m_num_points = 0;
            int m_num_lines = 0;
            int m_num_surfaces = 0;
            int m_num_volumes = 0;
            std::vector<point> m_points;
            std::map<int, int> m_point_map_index2tag;
            std::map<int, int> m_point_map_tag2index;
            std::vector<line> m_lines;
            std::map<int, int> m_line_map_index2tag;
            std::map<int, int> m_line_map_tag2index;
            std::vector<surface> m_surfaces;
            std::map<int, int> m_surface_map_index2tag;
            std::map<int, int> m_surface_map_tag2index;
            std::vector<volume> m_volumes;
            std::map<int, int> m_volume_map_index2tag;
            std::map<int, int> m_volume_map_tag2index;
    }; //geo
}

#endif // GMESH_H