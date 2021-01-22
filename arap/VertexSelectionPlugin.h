#ifndef IGL_OPENGL_GFLW_IMGUI_IMGUIDRAWLISTPLUGIN_H
#define IGL_OPENGL_GFLW_IMGUI_IMGUIDRAWLISTPLUGIN_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/AABB.h>
#include <igl/igl_inline.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <vector>

const unsigned char MOUSE_LEFT = 0;
const unsigned char MOUSE_MIDDLE = 1;
const unsigned char MOUSE_RIGHT = 2;

const int MODIFIER_SHIFT = 1;
const int MODIFIER_CTRL = 2;
const int MODIFIER_ALT = 4;

// 2 * MOUSE_AREA as a square is searched around a mouse click for a projected vertex
const int MOUSE_AREA = 8;

class Viewer;

class VertexSelectionPlugin : public igl::opengl::glfw::ViewerPlugin {
  public:
    // callback called when a vertex is clicked, vertex id is passed
    std::function<void(int)> callback_anchor_selected;
    std::function<void(int, Eigen::Vector3d)> callback_vertex_dragged;
    VertexSelectionPlugin() = default;
    void reset();

    IGL_INLINE bool mouse_down(int button, int modifier) override;
    IGL_INLINE bool mouse_up(int button, int modifier) override;
    IGL_INLINE bool mouse_move(int mouse_x, int mouse_y) override;

    std::set<int> fixedPoints;
  private:
    std::vector<Eigen::RowVector2f> L;
    Eigen::VectorXd W;
    Eigen::Array<double, Eigen::Dynamic, 1> and_visible;
    bool dragging = false;
    int drag_vertex_idx;

    int get_vertex_from_click(int mouse_x, int mouse_y);
};

#endif