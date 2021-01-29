#include "VertexSelectionPlugin.h"
#include <igl/screen_space_selection.h>

int VertexSelectionPlugin::get_vertex_from_click(int mouse_x, int mouse_y) {
    this->L.clear();
    this->L.emplace_back(mouse_x - MOUSE_AREA,
                         this->viewer->core().viewport(3) - mouse_y - MOUSE_AREA);
    this->L.emplace_back(mouse_x - MOUSE_AREA,
                         this->viewer->core().viewport(3) - mouse_y + MOUSE_AREA);
    this->L.emplace_back(mouse_x + MOUSE_AREA,
                         this->viewer->core().viewport(3) - mouse_y + MOUSE_AREA);
    this->L.emplace_back(mouse_x + MOUSE_AREA,
                         this->viewer->core().viewport(3) - mouse_y - MOUSE_AREA);

    this->W = Eigen::VectorXd::Zero(this->viewer->data().V.rows());
    this->and_visible = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(this->viewer->data().V.rows());

    igl::AABB<Eigen::MatrixXd, 3> tree;
    tree.init(this->viewer->data().V, this->viewer->data().F);

    screen_space_selection(this->viewer->data().V, this->viewer->data().F, tree, this->viewer->core().view,
                           this->viewer->core().proj, this->viewer->core().viewport, L, W, and_visible);

    for (int i = 0; i < W.size(); i++) {
        if (abs(W[i]) >= 0.5 && abs(and_visible[i]) >= 0.5) {
            return i;
        }
    }
    return -1;
}

void VertexSelectionPlugin::reset() {
    this->fixedPoints.clear();
    this->dragging = false;
}

bool VertexSelectionPlugin::mouse_down(int button, int modifier) {
    if (button != MOUSE_LEFT) {
        return false;
    }

    if ((modifier & MODIFIER_CTRL) == MODIFIER_CTRL) {
        auto selected = get_vertex_from_click(this->viewer->current_mouse_x, this->viewer->current_mouse_y);
        if (selected >= 0) {
            this->fixedPoints[selected] = this->viewer->data().V.row(selected);

            this->callback_anchor_selected(selected, this->viewer->data().V.row(selected));
        }
        return true;
    } else if ((modifier & MODIFIER_SHIFT) == MODIFIER_SHIFT) {
        auto selected = get_vertex_from_click(this->viewer->current_mouse_x, this->viewer->current_mouse_y);
        if (selected >= 0) {
            this->dragging = true;
            this->drag_vertex_idx = selected;
            this->fixedPoints[selected] = this->viewer->data().V.row(selected);
        }
        return true;
    }
    return false;
}

bool VertexSelectionPlugin::mouse_move(int button, int modifier) {
    //if (button == MOUSE_LEFT && (modifier & MODIFIER_SHIFT) == MODIFIER_SHIFT && this->dragging) {
    if (this->dragging) {
        //std::cout << "done some dragging" << std::endl;
        //this->dragging = false;

        Eigen::Vector4d vertex_homog = Eigen::Vector4d::Ones();
        vertex_homog.segment<3>(0) = this->viewer->data().V.row(this->drag_vertex_idx);
        // transform the selected vertex into camera coordinates
        Eigen::Vector4d vertex_in_camera = this->viewer->core().view.cast<double>() * vertex_homog;
        Eigen::Vector4d screen_pos = Eigen::Vector4d::Ones();
        screen_pos(0) = this->viewer->current_mouse_x / (this->viewer->core().viewport(2) / 2) - 1.0;
        screen_pos(1) = (this->viewer->core().viewport(3) - this->viewer->current_mouse_y) / (this->viewer->core().viewport(3) / 2) - 1.0;
        screen_pos = this->viewer->core().proj.inverse().cast<double>() * screen_pos;
        // scale by the z-position of the original vertex
        screen_pos *= -vertex_in_camera(2);
        // transform back to world coordinates
        screen_pos(3) = 1.0;
        screen_pos = this->viewer->core().view.inverse().cast<double>() * screen_pos;

        this->callback_vertex_dragged(this->drag_vertex_idx, screen_pos.segment<3>(0));
        this->fixedPoints[this->drag_vertex_idx] = screen_pos.segment<3>(0);

        return true;
    }
    return false;
}

bool VertexSelectionPlugin::mouse_up (int button, int modifier) {
    if (button == MOUSE_LEFT) {
        if (this->dragging) {
            this->callback_anchor_selected(this->drag_vertex_idx, this->viewer->data().V.row(this->drag_vertex_idx));
        }
        this->dragging = false;
    }
    return false;
}

/*
bool VertexSelectionPlugin::mouse_move(int mouse_x, int mouse_y) {
    if (this->dragging) {
        return true;
    }
    return false;
}*/
