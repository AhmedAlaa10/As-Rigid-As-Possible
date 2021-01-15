#include "Compute.h"

#include <vector>
#include <iostream>

#include "igl/slice.h"
#include <igl/polar_svd3x3.h>

Compute::Compute(const Eigen::MatrixXd& vertices, const Eigen::VectorXi& fixedVertices,
	const Eigen::MatrixXi& faces, int maxIterations)
	: vertices_(vertices),
	fixedVertices_Index(fixedVertices),
	faces_(faces),
	maxIterations_(maxIterations)
{
	//compute # of free vertices
	int nFree = vertices_.rows() - fixedVertices_Index.size();

	//resize freeVertices_index
	freeVertices_Index.resize(nFree);

	//Fill in the vector of indices for free vertices
	int i = 0;
	int j = 0;

	for (int k = 0; k < vertices_.rows(); k++)
	{
		if (i < fixedVertices_.size() && k == fixedVertices_Index(i))
		{
			i++;
		}
		else
		{
			freeVertices_Index(j) = k;
				j++;
		}
	}
	/*
	*TODO: Check that the sizes of the fixedVertices_Index and freeVertices_index are correct
	*************************************************************************************
	*************************************************************************************
	***                                                                               ***
	***                                                                               ***
	***                                                                               ***
	***                                                                               ***
	***                                                                               ***
	*************************************************************************************
	*/

}

void Compute::ComputeARAP(const Eigen::MatrixXd& fixedVertices)
{
	//Using the alternating minimization strategy.
	//
	//Preprocess(fixedVertices);

	int iteration = 0;

	while (iteration < maxIterations_)
	{
		SingleIteration();
		iteration++;
	}
}

void Compute::ComputeRotations()
{
	//check equations (5) and (6) in the paper for computing rotations

	//total number of vertices
	int nVertices = vertices_.rows();

	//a vector of matrices to hold the products in equation (5) in the paper.
	//S = per-edge weights * rest edge * deformed edge, summed over all vertices.
	std::vector<Eigen::MatrixXd> S_matrix;

	for (int v_1 = 0; v_1 < nVertices; v_1++)
	{
		for (auto neighbor : neighbors_[v_1])
		{
			int v_2 = neighbor.first;

			Eigen::Vector3d restEdge = vertices_.row(v_1) - vertices_.row(v_2);
			Eigen::Vector3d deformedEdge = updatedVertices_.row(v_1) - updatedVertices_.row(v_2);

			S_matrix[v_1] += weight_.coeff(v_1, v_2) * restEdge * deformedEdge.transpose();
		}
	}

	for (int v = 0; v < nVertices; v++)
	{
		Eigen::Matrix3d rotation;

		//polar_svd3x3 computes the polar decomposition (U, V) of a matrix using SVD
		//recall equation (6):			R[i] = V[i] * U[i].transpose()
		//however, polar_svd3x3 outputs R[i] = U[i] * V[i].transpose()
		//therefore we take the transpose of the output rotation.

		igl::polar_svd3x3(S_matrix[v], rotation);
		rotations[v] = rotation.transpose();
	}

}

void Compute::ComputeVertices()
{
	//step 1: write down the fixed vertices in the updatedVertices_ matrix

	int nFixed = fixedVertices_Index.size();
	int nFree = freeVertices_Index.size();

	for (int v = 0; v < nFixed; v++)
	{
		updatedVertices_.row(fixedVertices_Index(v)) = fixedVertices_.row(v);
	}

	//step 2: Factorize the sparse Laplace-Beltrami operator

	//TODO: Write a member function to do this task, since we used it in different places
	//neither Eigen::sparseLU solver nor Eigen::sparseMatrix L_operator are defined here but in the child class ArapCompute
	//so we need to define it here

	Eigen::SparseMatrix<double> L_operator;
	igl::slice(weight_, freeVertices_Index, freeVertices_Index, L_operator);
	//L_operator *= -1.0;
	L_operator.makeCompressed();

	Eigen::SparseLU<Eigen::SparseMatrix<double>> sparse_solver;
	sparse_solver.compute(L_operator);

	//step 3: Define the R.H.S of equation(9)



}

Eigen::MatrixXd& Compute::ComputeVertixGradients() const
{

}

double Compute::ComputeVertixGradients_(int vertix, int dimension) const
{

}


















