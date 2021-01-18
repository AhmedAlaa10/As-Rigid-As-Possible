#include <igl/opengl/glfw/Viewer.h>
#include <igl/project.h>
#include <igl/readOFF.h>
#include <igl/winding_number.h>

#include "VertexSelectionPlugin.h"
#include "ARAP_Compute.h"
#include "Compute.h"

int main(int argc, char** argv) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    int currentVertexID;
    // inputs for the ARAP
    int maxIter;
    Eigen::Vector3d oldPosition;
    Eigen::Vector3d newPosition;


    // optional command line argument for different files (or paths)
    std::string path = argc == 2 ? argv[1] : "../data/bunny.off";
    
    // load the mesh in off format, and abort on any read errors (unfortunately, libigl does not seem to
    // provide any further hints as to what exactly went wrong, so this just prints out a generic usage hint)
    if (!igl::readOFF(path, V, F)) {
        std::cout << "usage: " << argv[0] << " [meshfile.off]" << std::endl;
        return 1;
    }

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
