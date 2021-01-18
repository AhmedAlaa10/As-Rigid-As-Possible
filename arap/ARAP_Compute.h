#ifndef _ARAP_COMPUTATION_H_
#define _ARAP_COMPUTATION_H_


#include "Eigen/SparseLU"
#include "Compute.h"

class ArapCompute {

public:
    // Class Constructor
	ArapCompute(const Eigen::MatrixXd& vertices, const Eigen::VectorXi& fixedVertices,
		const Eigen::MatrixXi& faces, const int maxIterations);


    // Function to compute the weights according to the paper
	//($w_{ij}$) uses the faces and the vertices

    void ComputeWeights();

    // Initialization for the starting guess of the vertex shift
    void NaiveLaplaceEditing();

	double ComputeEnergy() const;

	Eigen::Vector3d ComputeCotangent(int face_index) const;

	//in equation(7) in the paper, there are two weights basically,
	//the per-edge weight we've been dealing with all the time,
	//and the per-cell weight, which can be assigned as the Voronoi area of the cell
	//However, this area term cancels out and we're left with only the per-edge weights(i.e. cotangents)
	std::vector<double> ComputeFaceArea(const Eigen::MatrixXi& faces) const;

	void ComputeRightHandSide();

	void alternatingOptimization();

private:

    // helper function to map from the vertex to the edge
    const int VertexToEdge_map[3][2] = { {0, 1}, {1, 2}, {2, 0} };

    // helper function to compute the angle between to vectors
    const double angleBetweenVectors(const Eigen::Vector3d& a, const Eigen::Vector3d& b);

    // Computes the indices of all the neighbouring vectors for a given vector ID
    // and stores them in a vetor of size N(i).
    Eigen::VectorXi computeNeighbourVertices(int vertexID);

    // Computes the Laplace Beltrami Operator used e.g. in equation (9)
    void computeLaplaceBeltramiOperator();

    // Applies the most recent Vertex Rotation to each Cell
    void ComputeRotations();


	//The Laplace-Beltrami operator
	Eigen::SparseMatrix<double> L_operator;

	//Sparse Linear Solver
	//About which sparse solver to use check https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
	//sparseLU is better for problems (large or small) with irregular patterns
	Eigen::SparseLU<Eigen::SparseMatrix<double>> sparse_solver;

		//a total no. of vertices x 3 matrix
	//each row represents a vertix position
	const Eigen::MatrixXd vertices_;

	//stores the output( vertices positions ) of the computeARAP function
	Eigen::MatrixXd updatedVertices_;

	//used to cache vertices passed in computeARAP or Preprocess
	Eigen::MatrixXd fixedVertices_;

	//a no. of faces x 3 matrix
	//the ith row stores the indices of the vertices of the ith face
	const Eigen::MatrixXi faces_;

	//fixedVertices_Index is a vector of the indices of the vertices we want to fix during the deformation.
	//the fixed vertices are defined by the interacting user.
	const Eigen::VectorXi fixedVertices_Index;

	//freeVertices_Index is a vector of the indices of the free vertices.
	//freeVertices_Index = total no. of vertices - no. of fixed indices
	Eigen::VectorXi freeVertices_Index;

	//maximum no. of iterations used to solve the ARAP problem
	const int maxIterations_;

	//stores all neighborhoods of each vertex
	std::vector<vertixNeighbors> neighbors_;

	//A sparse matrix used to store weights
	//This is a nVertices x nVertices matrix with zero diagonal elements
	//since the element i-j s.t. i=j is a not an edge but a point
	//the element i-j s.t. i!=j is the edge between vertix(i) and vertix(j)
	Eigen::SparseMatrix<double> weight_;

	//store rotations for all vertices in a vector of matrices
	//This can be viewed as the matrix R = 3nx3
	std::vector<Eigen::Matrix3d> rotations;

	//store
	Eigen::MatrixXd RHS;
};

#endif // _ARAP_COMPUTATION_H_


#pragma once
