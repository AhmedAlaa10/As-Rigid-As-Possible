#ifndef _ARAP_COMPUTATION_H_
#define _ARAP_COMPUTATION_H_


#include "Eigen/SparseLU"
#include "Compute.h"
#include <unordered_map>

typedef std::unordered_map<int, int> Map;

class ArapCompute 
{

public:
    // Class Constructor
	ArapCompute( const Eigen::MatrixXd& vertices,
				 const Eigen::VectorXi& fixedVertices,
				 const Eigen::MatrixXi& faces,
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

	void alternatingOptimization();

	//Computes the energy gives a set of updated vertices positions and rotations.
	double ComputeEnergy();

	//A getter function: returns the updatedVertices_ matrix.
	const Eigen::MatrixXd &getUpdatedVertices() const;

private:

    //Computes the cotangents of the edges in the given face
    Eigen::Vector3d ComputeCotangent(int face_index) const;

    //Helper function to map from the vertex to the edge
    const int VertexToEdge_map[3][2] = { {0, 1}, {1, 2}, {2, 0} };

    //Helper function to compute the angle between to vectors
    double angleBetweenVectors(const Eigen::Vector3d& a, const Eigen::Vector3d& b) const;

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
	*  stores the new calculated vertices
	*  used in ComputeEnergy
	*  if the computed energy is smaller than the previous energy calculated,
	*  then updatedVertices_ = cachedVertices_
	*/
	Eigen::MatrixXd cacheVertices_;

	Eigen::MatrixXd fixedVertices_;

	/*
	*  a no. of faces x 3 matrix
	*  the ith row stores the indices of the vertices of the ith face
	*/
	const Eigen::MatrixXi faces_;

	/* 
	*  ixedVertices_Index is a vector of the indices of the vertices we want to fix during the deformation.
	*  the fixed vertices are defined by the interacting user.
	*/
	const Eigen::VectorXi fixedVertices_Index;

	/*
	*  freeVertices_Index is a vector of the indices of the free vertices.
	*  freeVertices_Index = total no. of vertices - no. of fixed indices
	*/
	Eigen::VectorXi freeVertices_Index;

	//Maximum no. of iterations used to solve the ARAP problem
	int maxIterations_;

	//Stores all neighborhoods of each vertex
	std::vector<std::vector<int>> NeighborList;

	//A vector of map to store neighbors
	std::vector<Map> Neighbors_;

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

inline const Eigen::MatrixXd& ArapCompute::getUpdatedVertices() const
{
    return updatedVertices_;
}

#endif // _ARAP_COMPUTATION_H_


#pragma once
