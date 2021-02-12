#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <chrono>

#include <igl/arap.h>

#include "VertexSelectionPlugin.h"
#include "ARAP_Compute.h"


enum ARAP_MODE {
    OWN,
    LIBIGL
};



int main(int argc, char** argv) {
    Eigen::MatrixXd V;
    Eigen::MatrixXd deformedV;
    Eigen::MatrixXi F;
    

    // inputs for the ARAP
    int maxIter = 20;
    
    std::string path = "../data/bunny.off";
    
    ARAP_MODE mode = OWN;
    
    // reading in arguments. I really miss Rust's serde library right now ...
    if (argc >= 2) {
        if (!strcmp(argv[1],"libigl")) mode = LIBIGL;
        else if (!strcmp(argv[1],"own")) mode = OWN;
        else if (argc == 2) path = argv[1];
        else if (argc != 3) {
            std::cout << "wrong usage!" << std::endl;
            return 1;
        }
        if (argc == 3) path = argv[2];
    }
    
    // optional command line argument for different files (or paths)
    //std::string path = argc == 2 ? argv[1] : "../data/bunny.off";
    
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
    
    
    ArapCompute ownarap(V, plugin.fixedPoints, F, 20);
    igl::ARAPData iglarap;
    
    
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
        
        switch (mode) {
            case OWN:
                // perform a full ARAP pass, display & set the current mesh to the result
                ownarap.set_fixpoints(plugin.fixedPoints);
                ownarap.alternatingOptimization(50);
                viewer.data().set_vertices(ownarap.getUpdatedVertices());
                deformedV = ownarap.getUpdatedVertices();
                break;
            case LIBIGL:
                i = 0;
                Eigen::VectorXi fixed (plugin.fixedPoints.size());
                for (auto iter = plugin.fixedPoints.begin(); iter != plugin.fixedPoints.end(); ++iter) {
                    fixed(i) = iter->first;
                    i++;
                }
                deformedV = V;
                igl::arap_solve (fixpointmatrix, iglarap, deformedV);
                viewer.data().set_vertices(deformedV);
                break;
        }
        
    };
    plugin.callback_vertex_drag_start = [&](int vertexID, const Eigen::Vector3d& new_position) {
        std::cout << "INFO: dragging vertex #" << vertexID << ", doing simple deformation" << std::endl;
        auto now = std::chrono::high_resolution_clock::now();
        
        switch (mode) {
            case OWN:
                // "preview" transformation. start with the current (already-deformed) mesh to make it look better
                ownarap.set_fixpoints(plugin.fixedPoints);
                ownarap.alternatingOptimization(0);
                viewer.data().set_vertices(ownarap.getUpdatedVertices());
                deformedV = ownarap.getUpdatedVertices();
                break;
            case LIBIGL:
                Eigen::VectorXi fixed (plugin.fixedPoints.size());
                Eigen::MatrixXd fixpointmatrix (plugin.fixedPoints.size(), 3);
                int i = 0;
                for (auto iter = plugin.fixedPoints.begin(); iter != plugin.fixedPoints.end(); ++iter) {
                    fixed(i) = iter->first;
                    fixpointmatrix.row(i) = iter->second.transpose();
                    i++;
                }
                igl::arap_precomputation (V,F,V.cols(),fixed, iglarap);
                    deformedV = V;
                igl::arap_solve (fixpointmatrix, iglarap, deformedV);
                viewer.data().set_vertices(deformedV);
                break;
        }
                
    };

    auto last_update = std::chrono::high_resolution_clock::now();
    // update function: currently dragging a vertex
    plugin.callback_vertex_dragged = [&](int vertexID, const Eigen::Vector3d& new_position) {
        std::cout << "INFO: dragging vertex #" << vertexID << ", doing simple deformation" << std::endl;
        auto now = std::chrono::high_resolution_clock::now();
        if (std::chrono::duration_cast<std::chrono::milliseconds>(now - last_update).count() > 0.01) { // 200ms
            std::cout << "elapsed time since last update: " << (now - last_update).count() << "s" << std::endl;
            switch (mode) {
                case OWN:
                    // "preview" transformation. start with the current (already-deformed) mesh to make it look better
                    ownarap.iterate();
                    ownarap.set_fixpoints(plugin.fixedPoints);
                    viewer.data().set_vertices(ownarap.getUpdatedVertices());
                    deformedV = ownarap.getUpdatedVertices();
                    break;
                case LIBIGL:
                    Eigen::VectorXi fixed (plugin.fixedPoints.size());
                    Eigen::MatrixXd fixpointmatrix (plugin.fixedPoints.size(), 3);
                    int i = 0;
                    for (auto iter = plugin.fixedPoints.begin(); iter != plugin.fixedPoints.end(); ++iter) {
                        fixed(i) = iter->first;
                        fixpointmatrix.row(i) = iter->second.transpose();
                        i++;
                    }
                    igl::arap_precomputation (V,F,V.cols(),fixed, iglarap);
                    deformedV = V;
                    igl::arap_solve (fixpointmatrix, iglarap, deformedV);
                    viewer.data().set_vertices(deformedV);
                    break;
            }
            last_update = now;
        }
    };

    // user pressed a key
    viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer& viewer, unsigned int key, int modifiers) -> bool {
        if (key == 'i' || key == 'I') {
            switch (mode) {
                case OWN:
                    std::cout << "computing arap now" << std::endl;
                    ownarap.iterate();
                    viewer.data().set_vertices(ownarap.getUpdatedVertices());
                    deformedV = ownarap.getUpdatedVertices();
                    break;
                case LIBIGL:
                    std::cout << "iterating not supported in libigl mode!" << std::endl;
                    break;
            }
            return true;
        } else if (key == 'r' || key == 'R') {
            // redo the initial size hack
            plugin.reset();
            if (mode == OWN) {
                ownarap.set_fixpoints(plugin.fixedPoints);
                ownarap.alternatingOptimization(1);
                deformedV = ownarap.getUpdatedVertices();
            } else {
                deformedV = V;
            }
            viewer.data().set_vertices(deformedV);
            viewer.data().clear_points();
        }

        return false;
    };


    switch (mode) {
        case OWN:
            // initial ARAP pass, to scale the mesh down
            ownarap.alternatingOptimization();
            // copy the initial mesh, so we can restore it later
            deformedV = ownarap.getUpdatedVertices();
            break;
        case LIBIGL:
            iglarap.energy = igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS; 
            deformedV = V;
            break;
    }
    
    viewer.plugins.push_back(&plugin);

    
    
    viewer.data().set_mesh(deformedV, F);
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    
    std::cout << "\n ==== ARAP keybindings ====\n"
              << "  Ctrl    Make fixpoint\n"
              << "  Shift   Drag fixpoint with mouse\n"
              << "  i       Iterate ARAP one step further\n"
              << "  r       Reset mesh" << std::endl;
    viewer.launch();
    
}
