#pragma once

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

};
/*
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
*/
