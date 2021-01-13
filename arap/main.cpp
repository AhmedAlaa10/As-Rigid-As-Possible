#include <igl/opengl/glfw/Viewer.h>
#include <igl/project.h>
#include <igl/readOFF.h>
#include <igl/winding_number.h>

#include "VertexSelectionPlugin.h"

int main() {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    // load the mesh in off format
    igl::readOFF("../data/bunny.off", V, F);

    // show the mesh in the igl viewer
    igl::opengl::glfw::Viewer viewer;
    VertexSelectionPlugin plugin;
    viewer.data().point_size = 20;

    plugin.callback_anchor_selected = [&](int vertexID) {
        std::cout << "selected an anchor point at " << V.row(vertexID) << std::endl;
        viewer.data().add_points(V.row(vertexID), Eigen::RowVector3d(255, 0, 0));
    };

    plugin.callback_vertex_dragged = [&](int vertexID, const Eigen::Vector3d& new_position) {
        std::cout << "updating vertex\n" << V.row(vertexID) << "\nto new vertex position\n" << new_position << std::endl;
        V.row(vertexID) = new_position;
        viewer.data().set_vertices(V);
    };

    viewer.plugins.push_back(&plugin);
    viewer.data().set_mesh(V, F);
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.launch();
}
