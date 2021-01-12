#include "VertexSelectionPlugin.h"
#include <igl/screen_space_selection.h>

bool VertexSelectionPlugin::mouse_down(int button, int modifier) {
    if (button == MOUSE_LEFT && (modifier & MODIFIER_CTRL) == MODIFIER_CTRL) {
        std::vector<Eigen::RowVector2f> L;
        L.emplace_back(this->viewer->current_mouse_x - MOUSE_AREA,
                       this->viewer->core().viewport(3) - this->viewer->current_mouse_y - MOUSE_AREA);
        L.emplace_back(this->viewer->current_mouse_x - MOUSE_AREA,
                       this->viewer->core().viewport(3) - this->viewer->current_mouse_y + MOUSE_AREA);
        L.emplace_back(this->viewer->current_mouse_x + MOUSE_AREA,
                       this->viewer->core().viewport(3) - this->viewer->current_mouse_y + MOUSE_AREA);
        L.emplace_back(this->viewer->current_mouse_x + MOUSE_AREA,
                       this->viewer->core().viewport(3) - this->viewer->current_mouse_y - MOUSE_AREA);

        Eigen::VectorXd W = Eigen::VectorXd::Zero(this->viewer->data().V.rows());
        Eigen::Array<double, Eigen::Dynamic, 1> and_visible =
                Eigen::Array<double, Eigen::Dynamic, 1>::Zero(this->viewer->data().V.rows());

        igl::AABB<Eigen::MatrixXd, 3> tree;
        tree.init(this->viewer->data().V, this->viewer->data().F);

        screen_space_selection(this->viewer->data().V, this->viewer->data().F, tree, this->viewer->core().view,
                               this->viewer->core().proj, this->viewer->core().viewport, L, W, and_visible);

        for (int i = 0; i < W.size(); i++) {
            if (abs(W[i]) >= 0.5 && abs(and_visible[i]) >= 0.5) {
                this->callback(i);
            }
        }

        return true;
    }
    return false;
}

bool VertexSelectionPlugin::mouse_up(int button, int modifier) {
    if (button == MOUSE_LEFT && (modifier & MODIFIER_CTRL) == MODIFIER_CTRL) {
        return true;
    }
    return false;
}

bool VertexSelectionPlugin::mouse_move(int mouse_x, int mouse_y) {
    return false;
}

