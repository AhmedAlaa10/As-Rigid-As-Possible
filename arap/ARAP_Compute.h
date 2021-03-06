#pragma once

#include "Eigen/SparseLU"
#include "Compute.h"
#include <unordered_map>
#include <set>

typedef std::unordered_map<int, int> Map;

class ArapCompute {

public:
    // Class Constructor
    ArapCompute(const Eigen::MatrixXd &vertices,
                const std::map<int, Eigen::Vector3d> &fixedVertices,
                const Eigen::MatrixXi &faces,
                const int maxIterations);

    /*
	* Function to compute the weights according to the paper
	* W_{ij}) uses the faces and the vertices
	*/
    void ComputeWeights();

    void UpdateVertices();

    //Computes the Laplace Beltrami Operator used e.g. in equation (9)
    void computeLaplaceBeltramiOperator();

    //Initialization for the starting guess of the vertex shift
    void NaiveLaplaceEditing();

    //Applies the most recent Vertex Rotation to each Cell
    void ComputeRotations();

    void ComputeRightHandSide();

    void alternatingOptimization(int maxIterations = 0);
    
    void set_fixpoints (const std::map<int, Eigen::Vector3d> &fixedVertices);
    
    void iterate();

    //A getter function: returns the updatedVertices_ matrix.
    const Eigen::MatrixXd &getUpdatedVertices() const;

private:

    //Computes the cotangents of the edges in the given face
    Eigen::Vector3d ComputeCotangent(int face_index) const;

    /*
	*  Computes the indices of all the neighbouring vectors for a given vector ID
    *  and stores them in a vetor of size N(i).
	*/
    void computeNeighbourVertices();

    //The Laplace-Beltrami operator
    Eigen::SparseMatrix<double> L_operator;

    /*
    *  Sparse Linear Solver
    *  About which sparse solver to use check https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
    *  sparseLU is better for problems (large or small) with irregular patterns
    */
    Eigen::SparseLU<Eigen::SparseMatrix<double>> sparse_solver;

    //A total no. of vertices x 3 matrix, each row represents a vertix position
    Eigen::MatrixXd vertices_;

    //Stores the output( vertices positions ) of the computeARAP function
    Eigen::MatrixXd updatedVertices_;

    /*
    *  a no. of faces x 3 matrix
    *  the ith row stores the indices of the vertices of the ith face
    */
    const Eigen::MatrixXi faces_;

    /*
    *  Fixed vertices as anchor points. Maps the vertex id to the new, fixed position
    */
    std::map<int, Eigen::Vector3d> fixedVertices;

    //Maximum no. of iterations used to solve the ARAP problem
    int maxIterations_;

    //Stores all neighborhoods of each vertex
    std::vector<std::set<int>> NeighborList;

    /*
    *  A sparse matrix used to store weights
    *  This is a nVertices x nVertices matrix with zero diagonal elements
    *  since the element i-j s.t. i=j is a not an edge but a point
    *  the element i-j s.t. i!=j is the edge between vertix(i) and vertix(j)
    */
    Eigen::SparseMatrix<double> weight_;

    /*
    *  store rotations for all vertices in a vector of matrices
    *  This can be viewed as the matrix R = 3nx3
    */
    std::vector<Eigen::Matrix3d> rotations;

    //Stores the right-hand-side of equation(9).
    Eigen::MatrixXd RHS;

};

inline const Eigen::MatrixXd &ArapCompute::getUpdatedVertices() const {
    return updatedVertices_;
}
