
#ifndef _COMPUTATION_H_
#define _COMPUTATION_H_

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include <unordered_map>

typedef std::unordered_map<int, int> vertixNeighbors;

class Compute {
public:

	Compute(const Eigen::MatrixXd& vertices, const Eigen::VectorXi& fixedVertices,
		const Eigen::MatrixXi& faces, int maxIterations);

	virtual ~Compute() {}

	//pre-factorizes per-edge weights (i.e. cotangents) and the left matrix
	virtual void Prefactorize();

	//this function does pre-processing before solving ARAP
	virtual void Preprocess(const Eigen::MatrixXd& fixed_vertices) = 0;

	//solves ARAP for a single iteration only
	//helps us analyze the algorithm
	virtual void SingleIteration();

	//computes the ARAP problem
	virtual void ComputeARAP(const Eigen::MatrixXd& fixedVertices);

	//returns updated vertices
	const Eigen::MatrixXd& GetNewVertices() const;

	//return fixed vertices
	const Eigen::VectorXi& GetFixedVertices() const;

	//returns faces
	const Eigen::MatrixXi& GetFaces() const;

	//returns maximum number of iterations
	const int GetIterations() const;

	//computes Energy
	const double ComputeEnergy() const;

	//computes rotations based on given vertices, as of equation (6) in the arap paper
	void ComputeRotations();

	//computes vertices based on given rotations
	void ComputeVertices();

	//computes the gradients of the energy w.r.t all vertix positions
	//returns a #no. of vertices x 3 matrix of gradients of vertices position w.r.t all dimensions
	Eigen::MatrixXd& ComputeVertixGradients() const;


	//Given a vertix position p' and the dimension (i.e. x, y or z)
	//returns the energy gradient w.r.t to it
	double ComputeVertixGradients_(int vertix, int dimension) const;






protected:


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

};

inline const Eigen::MatrixXd& Compute::GetNewVertices() const {
	return updatedVertices_;
}

inline const Eigen::MatrixXi& Compute::GetFaces() const {
	return faces_;
}

inline const Eigen::VectorXi& Compute::GetFixedVertices() const {
	return fixedVertices_Index;
}

inline const int Compute::GetIterations() const {
	return maxIterations_;
}


#endif














#pragma once
