#[[
# Choose which chapters to compile
option(LIBIGL_TUTORIALS_CHAPTER6 "Compile libigl tutorials chapter 6" ON)

# Download tutorial data
include(libigl_tutorial_data)

if(NOT TARGET igl::glfw)
    message(WARNING "Most tutorial executables depend on OpenGL and glfw. Use `cmake -DLIBIGL_GLFW=ON ..`")
endif()

# Chapter 6
if(LIBIGL_TUTORIALS_CHAPTER6)
    igl_add_tutorial(610_CSGTree igl::glfw igl_copyleft::cgal)
    igl_copy_dll(610_CSGTree)
endif()
]]

include(libigl_tutorial_data)
igl_add_tutorial(610_CSGTree igl::glfw igl_copyleft::cgal)
igl_copy_dll(610_CSGTree)
