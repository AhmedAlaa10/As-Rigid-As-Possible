#ifndef _ARAP_COMPUTATION_H_
#define _ARAP_COMPUTATION_H_

#include "Compute.h""

#include "Eigen/SparseLU"

class ArapCompute : public Compute {

public:

	ArapCompute(const Eigen::MatrixXd& vertices, const Eigen::VectorXi& fixedVertices,
		const Eigen::MatrixXi& faces, const int maxIterations);


	void Prefactorize();

	void Preprocess(const Eigen::MatrixXd& fixedVertices);

	void SingleIteration();

	void computeARAP();

	double ComputeEnergy() const;

	Eigen::Vector3d ComputeCotangent(int face_index) const;

	//in equation(7) in the paper, there are two weights basically,
	//the per-edge weight we've been dealing with all the time,
	//and the per-cell weight, which can be assigned as the Voronoi area of the cell
	//However, this area term cancels out and we're left with only the per-edge weights(i.e. cotangents)
	std::vector<double> ComputeFaceArea(const Eigen::MatrixXi& faces) const;

private:

	//The Laplace-Beltrami operator
	Eigen::SparseMatrix<double> L_operator;

	//Sparse Linear Solver
	//About which sparse solver to use check https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
	//sparseLU is better for problems (large or small) with irregular patterns
	Eigen::SparseLU<Eigen::SparseMatrix<double>> sparse_solver;
};

#endif _ARAP_COMPUTATION_H_


#pragma once
