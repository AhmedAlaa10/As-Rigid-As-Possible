#include "ARAP_Compute.h"
#include <iostream>

#include "igl/slice.h"
#include "igl/polar_svd.h"
#include "igl/polar_svd3x3.h"

using namespace std;


ArapCompute::ArapCompute(const Eigen::MatrixXd& vertices,
                         const Eigen::VectorXi& fixedVertices,
	                     const Eigen::MatrixXi& faces, 
                         int maxIterations)
	: vertices_(vertices),
      faces_(faces),
      fixedVertices_Index(fixedVertices),
      maxIterations_(maxIterations){

    //Number of free vertices.
    int nVertices = vertices_.rows();
    int nFixed    = fixedVertices_Index.size();
    int nFree     = nVertices - nFixed;

    //Initialize the vector holding the indices of the free vertices.
    freeVertices_Index.resize(nFree);

    int i = 0, j = 0;
    for (int k = 0; k < nVertices; k++)
    {
        if (i < nFixed && k == fixedVertices_Index(i))
        {
            i++;
        }
        else
        {
            freeVertices_Index(j) = k;
            j++;
        }

    }
    //Test these for loops
    if (i != nFixed || j != nFree)
    {
        std::cout<<"There is a mistake in dimensions of freeVertices and/or fixedVertices"<<std::endl;
    }

}

void ArapCompute::ComputeWeights()
{
    std::cout << "Computing weights ..." << std::endl;

    //compute weights
    int nVertices = vertices_.rows();
    int nFaces = faces_.rows();

    //Initial the weight matrix
    weight_.resize(nVertices, nVertices);
    weight_.setZero();

    // iterating through each face to work with every single triangle
    for (int face = 0; face < nFaces; face++)
    {
        //compute the cotangent vector of the face
        Eigen::Vector3d cotan = ComputeCotangent(face);

        //loop over the three vertices within the face
        for (int v = 0; v < 3; v++)
        {
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

Eigen::Vector3d ArapCompute::ComputeCotangent(int face_index) const
{
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
    //let's call the sides a, b, c
    double v3_squared = (v1 - v2).squaredNorm();
    double v1_squared = (v2 - v3).squaredNorm();
    double v2_squared = (v3 - v1).squaredNorm();

    cotangent(0) = (v3_squared + v1_squared - v2_squared) / (4 * area);
    cotangent(1) = (v2_squared + v1_squared - v3_squared) / (4 * area);
    cotangent(2) = (v3_squared + v2_squared - v1_squared) / (4 * area);

    return cotangent;
}

double ArapCompute::angleBetweenVectors(const Eigen::Vector3d& a, const Eigen::Vector3d& b) const
{
    double angle = 0.0;

    angle = std::atan2(a.cross(b).norm(), a.dot(b));

    return angle;
}

void ArapCompute::computeNeighbourVertices()
{
    std::cout << "computing neighboring vertices" << std::endl;
    int nVertices = vertices_.rows();
    NeighborList.resize(nVertices);

    for (int face = 0; face < faces_.rows(); face++)
    {
        for (int i = 0; i < 3; i++)
        {
            int firstVertex  = faces_(face, VertexToEdge_map[i][0]);
            int secondVertex = faces_(face, VertexToEdge_map[i][1]);


            /*
            *  IMPORTANT
            *
            *  NeighborList is a vector of vector holding the neighbor of each vertex in the mesh
            *  For each vertex i there corresponds a vector of neighboring vertices.
            *  each edge in every face has two vertices right, say i and j.
            *  i and j are neighbor of one another. so i is a neighbor of j and vice versa.
            *  This is basically what those two line down there mean.            *
            *  For vector NeighborList at the row(i) we add the vertex (j) to it's vector of neighboring vertices.
            *
            *  Now, how do we call back the neighbor of vertex i when we need them? This is simple
            *  first iterate over the vertices: for(int i=0; i<number_of_vertices; i++)
            *  then iterate over the neighbors: for( int v : NeighborList[i])
            *  which gives every vertex v neighboring vertex i
            *
            *  I hope i you got this right.
            */
            NeighborList[firstVertex].push_back(secondVertex);
            NeighborList[secondVertex].push_back(firstVertex);
        }
    }

}

void ArapCompute::computeLaplaceBeltramiOperator()
{
    std::cout << "Compute Laplace-Beltrami Operator ..." << std::endl;

    int nFree = freeVertices_Index.size();
    double weight;

	//compute the laplace-beltrami operator
	L_operator.resize(nFree, nFree);

    //Iteratre over the free vertices
	for (int i=0; i < nFree; i++)
    {
        //get the index of the free vertex i
        int j = freeVertices_Index(i);

        //iterate over the neighbors of the vertix i
        for (auto& neighbor : NeighborList[j])
        {
            weight = weight_.coeff(j, neighbor);

            L_operator.coeffRef(i, i) += weight;
        }

	}

	//for reducing memory consumption
	//makeCompressed() turns the sparseMatrix L_operator into the Compressed row/col storage form.
	L_operator.makeCompressed();

    sparse_solver.compute(L_operator);

    if (sparse_solver.info() != Eigen::Success)
    {
        std::cout << "Solving the sparse Laplace-Beltrami Operator failed!" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "Laplace-Beltrami Operator computed ... DONE !" << std::endl;
}

void ArapCompute::NaiveLaplaceEditing() {
    // computes the first guess for the vertex transformation P'_0 that is the basis for the
    // initial rotation
    std::cout << "Compute Initial guess ..." << std::endl;

    int n = vertices_.rows();

    Eigen::SparseMatrix<double> I(n, n);
    I.setIdentity();

    auto L_trans_L = L_operator.transpose() * L_operator;
    sparse_solver.compute(L_trans_L);

    auto L_Operator_inv = sparse_solver.solve(I);

    Eigen::VectorXd delta = L_operator * vertices_;

    // solution to ||Lp'-delta||^2 : p' = (L.transpose * L).inverser()*L.transpose * delta
    updatedVertices_ = L_Operator_inv * L_operator.transpose() * delta;

    std::cout << "Naive Laplacian Editing ... DONE !" << std::endl;
}

void ArapCompute::ComputeRotations()
{
    std::cout << "Compute rotations ..." << std::endl;

	//check equations (5) and (6) in the paper for computing rotations

	//total number of vertices
	int nVertices = vertices_.rows();

    Eigen::VectorXi neighbours;
    Eigen::Vector3d restEdge;
    Eigen::Vector3d deformedEdge;

	//a vector of matrices to hold the products in equation (5) in the paper.
	//S = per-edge weights * rest edge * deformed edge, summed over all vertices.
	std::vector<Eigen::MatrixXd> S_matrix;

	for (int i = 0; i < nVertices; i++)
	{
        //Iterate over neighbors of vertex i
        for (int j : NeighborList[i])
        {
            restEdge     = vertices_.row(i) - vertices_.row(j);
            deformedEdge = updatedVertices_.row(i) - updatedVertices_.row(j);

            S_matrix[i] += weight_.coeff(i, j) * restEdge * deformedEdge.transpose();
        }

	}

	for (int v = 0; v < nVertices; v++)
	{
		Eigen::Matrix3d rotation;

		//polar_svd3x3 computes the polar decomposition (U, V) of a matrix using SVD
		//recall equation (6):			R[i] = V[i] * U[i].transpose()
		//however, polar_svd3x3 outputs R[i] = U[i] * V[i].transpose()
		//therefore we take the transpose of the output rotation.
        Eigen::Matrix3d v_Rotation = S_matrix[v];

		igl::polar_svd3x3(v_Rotation, rotation);
		rotations[v] = rotation.transpose();
	}
    std::cout << "Compute rotations ... DONE!" << std::endl;
}

void ArapCompute::ComputeRightHandSide()
{
    std::cout << "Compute the right hand side ..." << std::endl;

    Eigen::VectorXi neighbours;
    int nVertices = vertices_.rows();
    int nFree = freeVertices_Index.size();

    //Initialize the right hand side matrix of equations (8) and (9).
    RHS = Eigen::MatrixXd::Zero(nFree, 3);

    for (int i = 0; i < nFree; i++)
    {
        Eigen::Vector3d sum(0.0, 0.0, 0.0);

        //Iterate over neighbors of vertex i
        for (int j : NeighborList[i])
        {
            Eigen::Vector3d Pi_minus_Pj = vertices_.row(i) - vertices_.row(j);
            Eigen::Matrix3d Ri_plus_Rj = rotations[i] + rotations[j];
            double weight = weight_.coeff(i, j) / 2.0;

            sum += weight * (rotations[i] + rotations[j]) * (vertices_.row(i) - vertices_.row(j)).transpose();
        }

        RHS(i,0) = sum.x();
        RHS(i,1) = sum.y();
        RHS(i,2) = sum.z();
    }

    std::cout << "Compute right hand side ... DONE!" << std::endl;
}

double ArapCompute::ComputeEnergy()
{
    int nVertices = vertices_.rows();

    //Initialize the total energy
    double total_energy = 0.0;

    //Loop over all vertices
    for (int i = 0; i < nVertices; i++)
    {

        //Loop over all neighbors
        for (auto& neighbor: Neighbors_[i])
        {
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

void ArapCompute::UpdateVertices()
{
    //step 1: write down the fixed vertices in the updatedVertices_ matrix

    int nFixed = fixedVertices_Index.size();
    int nFree = freeVertices_Index.size();
    int nVertices = vertices_.rows();

    Eigen::SparseMatrix<double> I(nVertices, nVertices);
    I.setIdentity();

    auto L_operator_inv = sparse_solver.solve(I);

    updatedVertices_ = L_operator_inv * RHS;

    for (int v = 0; v < nFixed; v++)
    {
        updatedVertices_.row(fixedVertices_Index(v)) = fixedVertices_.row(v);
    }
}

void ArapCompute::alternatingOptimization() {
    std::cout << "Alternating optimization ..." << std::endl;

    ComputeWeights();
    computeNeighbourVertices();
    computeLaplaceBeltramiOperator();
    NaiveLaplaceEditing();

    for (int iter = 0; iter < maxIterations_; iter++)
    {
        ComputeRotations();
        ComputeRightHandSide();
        UpdateVertices();

        vertices_ = updatedVertices_;
    }

    std::cout << "Optimization ... DONE !" << std::endl;
}