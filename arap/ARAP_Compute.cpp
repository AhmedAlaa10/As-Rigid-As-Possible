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
//    for (const auto& vertex : fixedVertices) {
//        updatedVertices_.row(vertex.first) = vertex.second;
//    }
}

void ArapCompute::ComputeWeights() {
    std::cout << "Computing weights ..." << std::endl;

    //compute weights
    int nVertices = vertices_.rows();
    int nFaces = faces_.rows();

    //Initial the weight matrix
    weight_.resize(nVertices, nVertices);
    weight_.setZero();

    // iterating through each face to work with every single triangle
    for (int face = 0; face < nFaces; face++) {
        //compute the cotangent vector of the face
        Eigen::Vector3d cotan = ComputeCotangent(face);

        //loop over the three vertices within the face
        for (int v = 0; v < 3; v++) {
            //indices of the two vertices in the edge
            int firstVertix = faces_(face, VertexToEdge_map[v][0]);
            int secondVertix = faces_(face, VertexToEdge_map[v][1]);

            //in Eigen::SparseMatrix, coeffRef() returns a non-constant reference to the value
            //of the matrix at position i, j
//            weight_.coeffRef(firstVertix, secondVertix) += 1.0 / tan(cotan(v)) / 2.0;
            weight_.coeffRef(firstVertix, secondVertix) += cotan(v) / 2.0;

            //since the weight of the edge i-j (edge between vertix(i) and vertix(j) is the same
            //of the edge j-i, then we assign them the same weight in the weight_ matrix
//            weight_.coeffRef(secondVertix, firstVertix) += 1.0 / tan(cotan(v)) / 2.0;
        }
    }

    std::cout << "Computing weights ... DONE !" << std::endl;
}

Eigen::Vector3d ArapCompute::ComputeCotangent(int face_index) const {
    //inititalize the cotangent vector
    //REMINDER: the cotangent is the double the per-edge weight in the paper (i.e. Wij = 1/2 * cotangent(i)
    Eigen::Vector3d cotangent(0.0, 0.0, 0.0);

    //3D positions of the vertices
    Eigen::Vector3d v1 = vertices_.row(faces_(face_index, 0));
    Eigen::Vector3d v2 = vertices_.row(faces_(face_index, 1));
    Eigen::Vector3d v3 = vertices_.row(faces_(face_index, 2));

    //Area of the triangle
    double area = (v1 - v2).cross(v2 - v3).norm() / 2;

    //the length of the sides squared
    //let's call the e1, e2, e3
    double e3_squared = (v1 - v2).squaredNorm();
    double e1_squared = (v2 - v3).squaredNorm();
    double e2_squared = (v3 - v1).squaredNorm();

    cotangent(0) = (e3_squared + e2_squared - e1_squared) / (4 * area);
    cotangent(1) = (e3_squared + e1_squared - e2_squared) / (4 * area);
    cotangent(2) = (e2_squared + e1_squared - e3_squared) / (4 * area);

    return cotangent;
//    Eigen::Vector3d cotangent(0.0, 0.0, 0.0);
//
//    //3D positions of the vertices
//    Eigen::Vector3d v0 = vertices_.row(faces_(face_index, 0));
//    Eigen::Vector3d v1 = vertices_.row(faces_(face_index, 1));
//    Eigen::Vector3d v2 = vertices_.row(faces_(face_index, 2));
//
//    Eigen::Vector3d e01 = v0 - v1;
//    Eigen::Vector3d e12 = v1 - v2;
//    Eigen::Vector3d e02 = v0 - v2;
//    cotangent(0) = acos(e12.dot(e02) / (e12.norm() * e02.norm())); // angle opposite edge 0-1
//    cotangent(1) = acos(e01.dot(e02) / (e01.norm() * e02.norm())); // angle opposite edge 1-2
//    cotangent(2) = acos(e01.dot(e12) / (e01.norm() * e12.norm())); // angle opposite edge 0-2
//
//    return cotangent;
}

double ArapCompute::angleBetweenVectors(const Eigen::Vector3d &a, const Eigen::Vector3d &b) const {
    double angle = 0.0;

    angle = std::atan2(a.cross(b).norm(), a.dot(b));

    return angle;
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
    double weight;

    //compute the laplace-beltrami operator
    L_operator.resize(nVertices, nVertices);
    L_operator.setZero();

    // Iteratre over all vertices
    for (int i = 0; i < nVertices; i++) {
        //get the index of the free vertex i
        //iterate over the neighbors of the vertix i
        for (auto &neighbor : NeighborList[i]) {
//                weight = weight_.coeff(i, neighbor);
            weight = 1.0;

            L_operator.coeffRef(i, i) += weight;
            L_operator.insert(i, neighbor) = -weight;
        }
    }

    // adjust the system matrix
    for (const auto& fixedVertex : fixedVertices) {
        auto k = fixedVertex.first;

        // delete rows and columns from the system matrix
        for (int j = 0; j < nVertices; j++) {
            if (L_operator.coeff(k, j) != 0.0) {
                L_operator.coeffRef(k, j) = 0.0;
                L_operator.coeffRef(j, k) = 0.0;
            }
        }
        L_operator.coeffRef(k, k) = 1.0;
    }

    //for reducing memory consumption
    //makeCompressed() turns the sparseMatrix L_operator into the Compressed row/col storage form.
    L_operator.makeCompressed();

    std::cout << "Laplace-Beltrami Operator computed ... DONE !" << std::endl;
}

void ArapCompute::NaiveLaplaceEditing() {
    // computes the first guess for the vertex transformation P'_0 that is the basis for the
    // initial rotation
    std::cout << "Compute Initial guess ..." << std::endl;

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> svd(L_operator);

    Eigen::MatrixXd delta = L_operator * vertices_;
    for (const auto& fixedVertex : fixedVertices) {
        delta.row(fixedVertex.first) = fixedVertex.second;
    }

    updatedVertices_ = svd.solve(delta);
    if (svd.info() != Eigen::Success) {
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
//            auto weight = weight_.coeff(i, j);
            auto weight = 1.0;

            sum += weight * restEdge * deformedEdge.transpose();

        }
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(sum, Eigen::ComputeFullU | Eigen::ComputeFullV);

        Eigen::Matrix3d rotation = svd.matrixV() * svd.matrixU().transpose();
        if (rotation.determinant() == -1) {
              Eigen::Matrix3d diag = Eigen::Matrix3d::Identity();
              diag(2, 2) = -1;
              rotation = svd.matrixV() * diag * svd.matrixU().transpose();
        }
        //polar_svd3x3 computes the polar decomposition (U, V) of a matrix using SVD
        //recall equation (6):			R[i] = V[i] * U[i].transpose()
        //however, polar_svd3x3 outputs R[i] = U[i] * V[i].transpose()
        //therefore we take the transpose of the output rotation.
//        igl::polar_svd3x3(rotMatrix, rotation);
        rotations[i] = rotation;
//        if ((rotation - Eigen::Matrix3d::Identity()).norm() > std::numeric_limits<float>::epsilon()) {
//            std::cout << "rotation matrix not equal to identity: \n" << rotation << std::endl;
//        }
//        std::cout << "rotation i: " << i << "of value:\n" << rotation << std::endl;
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
        //Iterate over neighbors of vertex i
        for (int j : NeighborList[i]) {
//            double weight = weight_.coeff(i, j) / 2.0;
            double weight = 1.0;

            RHS.row(i) -=
                    weight / 2.0 * (rotations[i] + rotations[j]) * (vertices_.row(i) - vertices_.row(j)).transpose();
        }
    }

    for (const auto& fixedVertex : fixedVertices) {
        auto k = fixedVertex.first;

        for (int j = 0; j < nVertices; j++) {
            if (L_operator.coeff(j, k) != 0.0) {
                RHS.row(j) -= L_operator.coeff(j, k) * fixedVertex.second;
            }
        }
        RHS.row(k) = fixedVertex.second;
    }

    std::cout << "Compute right hand side ... DONE!" << std::endl;
}

void ArapCompute::UpdateVertices() {
    //step 1: write down the fixed vertices in the updatedVertices_ matrix
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> svd(L_operator);
    updatedVertices_ = svd.solve(RHS);
}

void ArapCompute::alternatingOptimization() {
    std::cout << "Alternating optimization ..." << std::endl;

    ComputeWeights();
    computeNeighbourVertices();
    computeLaplaceBeltramiOperator();
//    NaiveLaplaceEditing();

    for (int iter = 0; iter < maxIterations_; iter++) {
        ComputeRotations();
        ComputeRightHandSide();
        UpdateVertices();
    }

    std::cout << "Optimization ... DONE !" << std::endl;
}
