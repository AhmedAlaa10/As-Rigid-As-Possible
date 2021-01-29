#include "ARAP_Compute.h"
#include <iostream>

#include "igl/slice.h"
#include "igl/polar_svd.h"

using namespace std;


ArapCompute::ArapCompute(const Eigen::MatrixXd &vertices,
                         const std::map<int, Eigen::Vector3d> &fixedVertices,
                         const Eigen::MatrixXi &faces,
                         int maxIterations)
        : vertices_(vertices),
          faces_(faces),
          fixedVertices(fixedVertices),
          maxIterations_(maxIterations),
          rotations(vertices.rows()),
          updatedVertices_(vertices) {

    // Initialize the vector holding the indices of the free vertices.
    for (const auto &vertex : fixedVertices) {
        updatedVertices_.row(vertex.first) = vertex.second;
    }
}

void ArapCompute::ComputeWeights() {
    std::cout << "Computing weights ..." << std::endl;

    //compute weights
    int nVertices = vertices_.rows();
    int nFaces = faces_.rows();

    //Initial the weight matrix
    weight_.resize(nVertices, nVertices);
    weight_.setZero();
//    const int vertexToEdge[3][2] = {{1, 2},
//                                    {0, 2},
//                                    {0, 1}};
    const int vertexToEdge[3][2] = {{0, 1},
                                    {1, 2},
                                    {2, 0}};

    // iterating through each face to work with every single triangle
//    for (int face = 0; face < nFaces; face++) {
//        //compute the cotangent vector of the face
//        //loop over the three vertices within the face
//        for (int v = 0; v < 3; v++) {
//            int v0 = faces_(face, v);
//            //indices of the two vertices in the edge
//            int v1 = faces_(face, vertexToEdge[v][0]);
//            int v2 = faces_(face, vertexToEdge[v][1]);
//
//            Eigen::Vector3d x = vertices_.row(v1) - vertices_.row(v0);
//            Eigen::Vector3d y = vertices_.row(v2) - vertices_.row(v0);
//            double w = x.dot(y) / (x.cross(y)).norm();
//            weight_.coeffRef(v1, v2) += w / 2.0;
//            weight_.coeffRef(v2, v1) += w / 2.0;
//        }
//
//    }
    for (int face = 0; face < nFaces; face++) {
        //compute the cotangent vector of the face
        Eigen::Vector3d cotan = ComputeCotangent(face);

        //loop over the three vertices within the face
        for (int v = 0; v < 3; v++) {
            //indices of the two vertices in the edge
            int firstVertix = faces_(face, vertexToEdge[v][0]);
            int secondVertix = faces_(face, vertexToEdge[v][1]);

            //in Eigen::SparseMatrix, coeffRef() returns a non-constant reference to the value
            //of the matrix at position i, j
            weight_.coeffRef(firstVertix, secondVertix) += cotan(v) / 2.0;
            weight_.coeffRef(secondVertix, firstVertix) += cotan(v) / 2.0;
        }
    }

    std::cout << "Computing weights ... DONE !" << std::endl;
}

Eigen::Vector3d ArapCompute::ComputeCotangent(int face_index) const {
    //inititalize the cotangent vector
    //REMINDER: the cotangent is the double the per-edge weight in the paper (i.e. Wij = 1/2 * cotangent(i)
//    Eigen::Vector3d cotangent(0.0, 0.0, 0.0);
//
//    //3D positions of the vertices
//    Eigen::Vector3d v1 = vertices_.row(faces_(face_index, 0));
//    Eigen::Vector3d v2 = vertices_.row(faces_(face_index, 1));
//    Eigen::Vector3d v3 = vertices_.row(faces_(face_index, 2));
//
//    //Area of the triangle
//    double area = (v1 - v2).cross(v2 - v3).norm() / 2;
//
//    //the length of the sides squared
//    //let's call the e1, e2, e3
//    double e3_squared = (v1 - v2).squaredNorm();
//    double e1_squared = (v2 - v3).squaredNorm();
//    double e2_squared = (v3 - v1).squaredNorm();
//
//    cotangent(0) = (e3_squared + e2_squared - e1_squared) / (4 * area);
//    cotangent(1) = (e3_squared + e1_squared - e2_squared) / (4 * area);
//    cotangent(2) = (e2_squared + e1_squared - e3_squared) / (4 * area);
//
//    return cotangent;
    Eigen::Vector3d cotangent(0.0, 0.0, 0.0);

    //3D positions of the vertices
    Eigen::Vector3d v0 = vertices_.row(faces_(face_index, 0));
    Eigen::Vector3d v1 = vertices_.row(faces_(face_index, 1));
    Eigen::Vector3d v2 = vertices_.row(faces_(face_index, 2));

    Eigen::Vector3d e01 = v0 - v1;
    Eigen::Vector3d e12 = v1 - v2;
    Eigen::Vector3d e02 = v0 - v2;
    cotangent(0) = acos(e12.dot(e02) / (e12.norm() * e02.norm())); // angle opposite edge 0-1
    cotangent(1) = acos(e01.dot(e02) / (e01.norm() * e02.norm())); // angle opposite edge 1-2
    cotangent(2) = acos(e01.dot(e12) / (e01.norm() * e12.norm())); // angle opposite edge 0-2

    return cotangent;
}

void ArapCompute::computeNeighbourVertices() {
    std::cout << "computing neighboring vertices" << std::endl;
    int nVertices = vertices_.rows();
    NeighborList.resize(nVertices);

    for (int face = 0; face < faces_.rows(); face++) {
        int v0 = faces_(face, 0);
        int v1 = faces_(face, 1);
        int v2 = faces_(face, 2);

        NeighborList[v0].emplace(v1);
        NeighborList[v0].emplace(v2);

        NeighborList[v1].emplace(v0);
        NeighborList[v1].emplace(v2);

        NeighborList[v2].emplace(v0);
        NeighborList[v2].emplace(v1);
    }
}

void ArapCompute::computeLaplaceBeltramiOperator() {
    std::cout << "Compute Laplace-Beltrami Operator ..." << std::endl;

    int nVertices = vertices_.rows();

    //compute the laplace-beltrami operator
    L_operator.resize(nVertices, nVertices);
    L_operator.setZero();

    // Iteratre over all vertices
    for (int i = 0; i < nVertices; i++) {
        if (fixedVertices.count(i) > 0) {
            L_operator.coeffRef(i, i) = 1.0;
        } else {
            //get the index of the free vertex i
            //iterate over the neighbors of the vertix i
            double diagWeight = 0;
            for (auto &neighbor : NeighborList[i]) {
                double weight = weight_.coeff(i, neighbor);
//                double weight = 1.0;

                diagWeight += weight;
                L_operator.insert(i, neighbor) = -weight;
            }
            L_operator.coeffRef(i, i) = diagWeight;
        }
    }

    //for reducing memory consumption
    //makeCompressed() turns the sparseMatrix L_operator into the Compressed row/col storage form.
    L_operator.makeCompressed();

    sparse_solver.analyzePattern(L_operator);
    sparse_solver.factorize(L_operator);

    std::cout << "Laplace-Beltrami Operator computed ... DONE !" << std::endl;
}

void ArapCompute::NaiveLaplaceEditing() {
    // computes the first guess for the vertex transformation P'_0 that is the basis for the
    // initial rotation
    std::cout << "Compute Initial guess ..." << std::endl;

    Eigen::MatrixXd delta = L_operator * vertices_;
    for (const auto &fixedVertex : fixedVertices) {
        delta.row(fixedVertex.first) = fixedVertex.second;
    }

    updatedVertices_ = sparse_solver.solve(delta);
    if (sparse_solver.info() != Eigen::Success) {
        std::cout << "Solving the sparse Laplace-Beltrami Operator failed!" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "Naive Laplacian Editing ... DONE !" << std::endl;
}

void ArapCompute::ComputeRotations() {
    std::cout << "Compute rotations ..." << std::endl;

    //check equations (5) and (6) in the paper for computing rotations

    //total number of vertices
    int nVertices = vertices_.rows();

    Eigen::VectorXi neighbours;
    Eigen::Vector3d restEdge;
    Eigen::Vector3d deformedEdge;

    //a vector of matrices to hold the products in equation (5) in the paper.
    //S = per-edge weights * rest edge * deformed edge, summed over all vertices.

    for (int i = 0; i < nVertices; i++) {
        Eigen::Matrix3d sum = Eigen::Matrix3d::Zero();
        //Iterate over neighbors of vertex i
        for (int j : NeighborList[i]) {
            restEdge = vertices_.row(i) - vertices_.row(j);
            deformedEdge = updatedVertices_.row(i) - updatedVertices_.row(j);
            auto weight = weight_.coeff(i, j);
//            auto weight = 1.0;

            sum += weight * restEdge * deformedEdge.transpose();

        }
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(sum, Eigen::ComputeFullU | Eigen::ComputeFullV);

        Eigen::Vector3d det(1, 1, (svd.matrixV() * svd.matrixU()).transpose().determinant());
        rotations[i] = svd.matrixV() * det.asDiagonal() * svd.matrixU().transpose();
        if (rotations[i].determinant() <= 0) {
            std::cout << "Error: determinant of rotation matrix <= 0" << std::endl;
        }
    }

    std::cout << "Compute rotations ... DONE!" << std::endl;
}

void ArapCompute::ComputeRightHandSide() {
    std::cout << "Compute the right hand side ..." << std::endl;

    Eigen::VectorXi neighbours;
    int nVertices = vertices_.rows();

    //Initialize the right hand side matrix of equations (8) and (9).
    RHS = Eigen::MatrixXd::Zero(nVertices, 3);
    for (int i = 0; i < nVertices; i++) {
        if (fixedVertices.count(i) > 0) {
            RHS.row(i) = fixedVertices.at(i);
        } else {
            //Iterate over neighbors of vertex i
            for (int j : NeighborList[i]) {
                double weight = weight_.coeff(i, j) / 2.0;
//                double weight = 1.0;

                RHS.row(i) +=
                        weight / 2.0 * (rotations[i] + rotations[j]) *
                        (vertices_.row(i) - vertices_.row(j)).transpose();
            }
        }
    }

    std::cout << "Compute right hand side ... DONE!" << std::endl;
}

void ArapCompute::UpdateVertices() {
    //step 1: write down the fixed vertices in the updatedVertices_ matrix
    updatedVertices_ = sparse_solver.solve(RHS);
}

void ArapCompute::alternatingOptimization() {
    std::cout << "Alternating optimization ..." << std::endl;

    ComputeWeights();
    computeNeighbourVertices();
    computeLaplaceBeltramiOperator();
    NaiveLaplaceEditing();

    for (int iter = 0; iter < maxIterations_; iter++) {
        ComputeRotations();
        ComputeRightHandSide();
        UpdateVertices();
    }

    std::cout << "Optimization ... DONE !" << std::endl;
}
