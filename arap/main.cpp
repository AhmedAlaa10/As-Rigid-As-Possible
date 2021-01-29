#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>

#include "VertexSelectionPlugin.h"
#include "ARAP_Compute.h"

int main(int argc, char** argv) {
    Eigen::MatrixXd V;
    Eigen::MatrixXd deformedV;
    Eigen::MatrixXi F;

    // inputs for the ARAP
    int maxIter = 20;

    // optional command line argument for different files (or paths)
    std::string path = argc == 2 ? argv[1] : "../data/bunny.off";
    
    // load the mesh in off format, and abort on any read errors (unfortunately, libigl does not seem to
    // provide any further hints as to what exactly went wrong, so this just prints out a generic usage hint)
    if (!igl::readOFF(path, V, F)) {
        std::cout << "usage: " << argv[0] << " [meshfile.off]" << std::endl;
        return 1;
    }

    // center mesh
    Eigen::Vector3d mean = V.colwise().mean();
    std::cout << "Centering mesh to its mean:\n" << mean << std::endl;
    for (int i = 0; i < V.rows(); i++) {
        V.row(i) -= mean;
    }

    // show the mesh in the igl viewer
    igl::opengl::glfw::Viewer viewer;
    VertexSelectionPlugin plugin;
    viewer.data().point_size = 20;
    
    // a new point was selected (i.e. click or dragging stopped)
    plugin.callback_anchor_selected = [&](int vertexID, const Eigen::Vector3d& position) {
        std::cout << "INFO: selected an anchor point at:\n" << position << std::endl;
        
        // the eigen viewer does not have a point_remove()-member function, so instead we
        // have to re-assemble the entire set of fixed point each time in a matrix format
        // which eigen understands
        Eigen::MatrixXd fixpointmatrix (plugin.fixedPoints.size(), 3);
        int i = 0;
        for (auto iter = plugin.fixedPoints.begin(); iter != plugin.fixedPoints.end(); ++iter){
            fixpointmatrix.row(i) = iter->second.transpose();
            i++;
        }
        viewer.data().set_points(fixpointmatrix, Eigen::RowVector3d(255,0,0));
        
        // perform a full ARAP pass, display & set the current mesh to the result
        std::cout << "INFO: computing full arap" << std::endl;
        ArapCompute arap(V, plugin.fixedPoints, F, maxIter);
        arap.alternatingOptimization();
        viewer.data().set_vertices(arap.getUpdatedVertices());
        deformedV = arap.getUpdatedVertices();
    };

    // update function: currently dragging a vertex
    plugin.callback_vertex_dragged = [&](int vertexID, const Eigen::Vector3d& new_position) {
        std::cout << "INFO: dragging vertex #" << vertexID << ", doing simple deformation" << std::endl;
        
        // "preview" transformation. start with the current (already-deformed) mesh to make it look better
        ArapCompute arap(deformedV, plugin.fixedPoints, F, 0);
        arap.alternatingOptimization();
        viewer.data().set_vertices(arap.getUpdatedVertices());
    };

    // user pressed a key
    viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer& viewer, unsigned int key, int modifiers) -> bool {
        /*if (key == 'i' || key == 'I') {
            std::cout << "computing arap now" << std::endl;
            ArapCompute arap(V, plugin.fixedPoints, F, maxIter);
            arap.alternatingOptimization();
            viewer.data().set_vertices(arap.getUpdatedVertices());
            deformedV = arap.getUpdatedVertices();
            return true;
        } else*/ if (key == 'r' || key == 'R') {
            // redo the initial size hack
            plugin.reset();
            ArapCompute arap(V, plugin.fixedPoints, F, 1);
            arap.alternatingOptimization();
            deformedV = arap.getUpdatedVertices();
            viewer.data().set_vertices(deformedV);
            viewer.data().clear_points();
        }

        return false;
    };

    viewer.plugins.push_back(&plugin);
    
    // initial ARAP pass, to scale the mesh down
    ArapCompute arap(V, plugin.fixedPoints, F, 1);
    // copy the initial mesh, so we can restore it later
    arap.alternatingOptimization();
    deformedV = arap.getUpdatedVertices();
    
    //deformedV = V;
    viewer.data().set_mesh(deformedV, F);
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    
    std::cout << "\n ==== ARAP keybindings ====\n"
              << "  Ctrl    Make fixpoint\n"
              << "  Shift   Drag fixpoint with mouse\n"
              /*<< "  i       Apply ARAP algorithm\n"*/
              << "  r       Reset mesh" << std::endl;
    viewer.launch();
}
