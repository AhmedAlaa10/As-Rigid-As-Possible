#include <igl/opengl/glfw/Viewer.h>
#include <igl/project.h>
#include <igl/readOFF.h>
#include <igl/winding_number.h>

#include "VertexSelectionPlugin.h"
#include "ARAP_Compute.h"
#include "Compute.h"

int main() {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    int currentVertexID;
    // inputs for the ARAP
    int maxIter;
    Eigen::Vector3d oldPosition;
    Eigen::Vector3d newPosition;


    // load the mesh in off format
    igl::readOFF("../data/bunny.off", V, F);

    // show the mesh in the igl viewer
    igl::opengl::glfw::Viewer viewer;
    VertexSelectionPlugin plugin;
    //Arap ArapPlugin;
    viewer.data().point_size = 20;

    plugin.callback_anchor_selected = [&](int vertexID) {
        std::cout << "selected an anchor point at " << V.row(vertexID) << std::endl;
        viewer.data().add_points(V.row(vertexID), Eigen::RowVector3d(255, 0, 0));
        currentVertexID = vertexID;
    };
    oldPosition = V.row(currentVertexID);

    plugin.callback_vertex_dragged = [&](int vertexID, const Eigen::Vector3d& new_position) {
        std::cout << "updating vertex\n" << V.row(vertexID) << "\nto new vertex position\n" << new_position << std::endl;
        V.row(vertexID) = new_position;
        viewer.data().set_vertices(V);
        //newPosition = new_position;
    };
/*
TODO's
define object of type Arap_compute
gothrough: compute weights (y->get the indexing straight) -> compute neighbours(y) -> compute L -> compute Initial Guess -> compute initial Rotation -> alternating optimization -> apply final rotation to cells
\detla = L * p
initial guess: min||L*p' = \delta||^2
initialize primary variables
prefactorize (compute cells vertices, weights, prefactor system matrix, initialRotation)
compute new positions(=p_i[->solve Lp' = b]) and update Rotation(=R_i=V_iU_i^T) iteratively
apply final rotation
*/



    viewer.plugins.push_back(&plugin);
    viewer.data().set_mesh(V, F);
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.launch();
}
