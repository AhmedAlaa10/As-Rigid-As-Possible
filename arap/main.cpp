#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>


int main() {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    // load the mesh in off format
    igl::readOFF("../data/bunny.off", V, F);

    // show the mesh in the igl viewer
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.launch();
}