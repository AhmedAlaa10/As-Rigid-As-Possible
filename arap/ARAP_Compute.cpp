#include "ARAP_Compute.h"
#include <iostream>

#include "igl/slice.h"
#include "igl/polar_svd.h"
#include "igl/polar_svd3x3.h"

#include <Eigen/Dense>

using namespace std;


ArapCompute::ArapCompute(const Eigen::MatrixXd &vertices,
                         const std::vector<int> &fixedVertices,
                         const Eigen::MatrixXi &faces,
                         int maxIterations)
        : vertices_(vertices),
          faces_(faces),
          fixedVertices_Index(fixedVertices),
          maxIterations_(maxIterations),
          rotations(vertices.rows()),
          updatedVertices_(vertices.rows(), vertices.cols()),
          fixedVertices_(fixedVertices.size(), 3) {

    //Initialize the vector holding the indices of the free vertices.
    updatedVertices_ = vertices_;

    for (int i = 0; i < fixedVertices.size(); i++) {
        fixedVertices_.row(i) = vertices_.row(fixedVertices_Index[i]);
        updatedVertices_.row(i) = vertices_.row(fixedVertices_Index[i]);
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
            weight_.coeffRef(firstVertix, secondVertix) += cotan(v) / 2.0;

            //since the weight of the edge i-j (edge between vertix(i) and vertix(j) is the same
            //of the edge j-i, then we assign them the same weight in the weight_ matrix
            weight_.coeffRef(secondVertix, firstVertix) += cotan(v) / 2.0;
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
        for (int i = 0; i < 3; i++) {
            int firstVertex = faces_(face, VertexToEdge_map[i][0]);
            int secondVertex = faces_(face, VertexToEdge_map[i][1]);

            /*
            *  IMPORTANT
            *
            *  NeighborList is a vector of vector holding the neighbor of each vertex in the mesh
            *  For each vertex i there corresponds a vector of neighboring vertices.
            *  each edge in every face has two vertices right, say i and j.
            *  i and j are neighbor of one another. so i is a neighbor of j and vice versa.
            *  This is basically what those two lines down there mean.            *
            *  For vector NeighborList at the row(i) we add the vertex (j) to it's vector of neighbors.
            *
            *  Now, how do we call back the neighbors of vertex i when we need them? This is simple
            *  first iterate over the vertices: for(int i=0; i<number_of_vertices; i++)
            *  then iterate over the neighbors: for( int v : NeighborList[i])
            *  which gives every vertex v neighboring vertex i
            *
            *  I hope i you got this right.
            */
            NeighborList[firstVertex].push_back(secondVertex);
            //NeighborList[secondVertex].push_back(firstVertex);
            // commented this out, since otherwise each neighbour is twice in each list
            // (as far as i can tell, the result is actually correct, and (up to order) identical to
            // what we would get if we have both these calls prefixed with a check if that vertex is
            // already in this particular list)
        }
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
<<<<<<< HEAD
        if (std::find(fixedVertices_Index.begin(), fixedVertices_Index.end(), i) != fixedVertices_Index.end()) 
        {
            L_operator.coeffRef(i, i) = 1.0;
        }
        else 
        {
            for (auto &neighbor : NeighborList[i]) 
            {
=======
        for (auto &neighbor : NeighborList[i]) {
>>>>>>> a3663cf64f915655cd9fc6f71b89c88d96c4e39e
                weight = weight_.coeff(i, neighbor);
//            weight = 1.0;

            L_operator.coeffRef(i, i) += weight;
            L_operator.coeffRef(i, neighbor) = -weight;
        }
    }

    // adjust the system matrix
    for (int i = 0; i < fixedVertices_Index.size(); i++) {
        auto k = fixedVertices_Index[i];

        // delete rows and columns from the system matrix
        for (int j = 0; j < nVertices; j++) {
            if (L_operator.coeff(k, j) != 0.0) {
                L_operator.coeffRef(k, j) = 0.0;
                L_operator.coeffRef(j, k) = 0.0;
            }
        }
        L_operator.coeffRef(k, k) = 1.0;
    }
    L_operator.makeCompressed();

    std::cout << "Laplace-Beltrami Operator computed ... DONE !" << std::endl;
}

void ArapCompute::NaiveLaplaceEditing() {
    // computes the first guess for the vertex transformation P'_0 that is the basis for the
    // initial rotation
    std::cout << "Compute Initial guess ..." << std::endl;

<<<<<<< HEAD
    int n = vertices_.rows();

    Eigen::SparseLU<Eigen::SparseMatrix<double>> naive_laplace_solver;

    auto L_trans_L = L_operator.transpose() * L_operator;
    naive_laplace_solver.compute(L_trans_L);

    Eigen::SparseMatrix<double> I(n, n);
    I.setIdentity();

    auto L_trans_L_inverse = naive_laplace_solver.solve(I);
    


    Eigen::MatrixXd delta = L_operator * vertices_;

    // solution to ||Lp'-delta||^2 : p' = (L.transpose * L).inverser()*L.transpose * delta
    updatedVertices_ = L_trans_L_inverse * L_operator.transpose() * delta;
=======
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> svd(L_operator);

    Eigen::MatrixXd delta = L_operator * vertices_;
    updatedVertices_ = svd.solve(delta);
    if (svd.info() != Eigen::Success) {
        std::cout << "Solving the sparse Laplace-Beltrami Operator failed!" << std::endl;
        exit(EXIT_FAILURE);
    }
>>>>>>> a3663cf64f915655cd9fc6f71b89c88d96c4e39e

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
        Eigen::Matrix3d rotation;
        Eigen::Matrix3d rotMatrix = Eigen::Matrix3d::Zero();
        //Iterate over neighbors of vertex i
        for (int j : NeighborList[i]) {
            restEdge = vertices_.row(i) - vertices_.row(j);
            deformedEdge = updatedVertices_.row(i) - updatedVertices_.row(j);
            auto weight = weight_.coeff(i, j);
//            auto weight = 1.0;

            rotMatrix += weight * restEdge * deformedEdge.transpose();

            //polar_svd3x3 computes the polar decomposition (U, V) of a matrix using SVD
            //recall equation (6):			R[i] = V[i] * U[i].transpose()
            //however, polar_svd3x3 outputs R[i] = U[i] * V[i].transpose()
            //therefore we take the transpose of the output rotation.
            igl::polar_svd3x3(rotMatrix, rotation);
            rotations[i] = rotation.transpose();
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
        //Iterate over neighbors of vertex i
        for (int j : NeighborList[i]) {
            double weight = weight_.coeff(i, j) / 2.0;
//            double weight = 1.0;

            RHS.row(i) +=
                    weight / 2.0f * (rotations[i] + rotations[j]) * (vertices_.row(i) - vertices_.row(j)).transpose();
        }
    }


    for (int i = 0; i < fixedVertices_Index.size(); i++) {
        auto k = fixedVertices_Index[i];
        auto fixedVertex = fixedVertices_.row(i);

        for (int j = 0; j < nVertices; j++) {
            if (L_operator.coeff(j, k) != 0.0) {
                RHS.row(i) -= L_operator.coeff(j, k) * fixedVertex;
            }
        }
        RHS.row(k) = fixedVertex;
    }

    std::cout << "Compute right hand side ... DONE!" << std::endl;
}

double ArapCompute::ComputeEnergy() {
    int nVertices = vertices_.rows();

    //Initialize the total energy
    double total_energy = 0.0;

    //Loop over all vertices
    for (int i = 0; i < nVertices; i++) {

        //Loop over all neighbors
        for (auto &neighbor: Neighbors_[i]) {
            int j = neighbor.first;

            double weight = weight_.coeff(i, j);
            double edge_energy = 0.0;

            //See equation (7) in the paper
            Eigen::Vector3d updatedVertices_vector = (cacheVertices_.row(i) - cacheVertices_.row(j)).transpose();
            Eigen::Vector3d oldVertices_vector = (vertices_.row(i) - vertices_.row(j)).transpose();

            total_energy += weight * (updatedVertices_vector - rotations[i] * oldVertices_vector).squaredNorm();
        }
    }
    return total_energy;
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
