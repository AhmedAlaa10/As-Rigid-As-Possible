#include "ARAP_Compute.h"
#include <iostream>

#include "igl/slice.h"
#include "igl/polar_svd.h"
#include "igl/polar_svd3x3.h"

using namespace std;


ArapCompute::ArapCompute(const Eigen::MatrixXd& vertices,
                         const Eigen::VectorXi& fixedVertices,
	                     const Eigen::MatrixXi& faces, 
                         const int maxIterations)
	: vertices_(vertices),
      faces_(faces),
      fixedVertices_(fixedVertices),
      maxIterations_(maxIterations){}

void ArapCompute::ComputeWeights()
{
    std::cout << "Computing weights ..." << std::endl;

    //compute weights
    int nVertices = vertices_.rows();
    int nFaces = faces_.size();

    int firstVertex;
    int secondVertex;
    int thirdVertex;

    Eigen::Vector3d A;
    Eigen::Vector3d B;
    Eigen::Vector3d C;

	weight_.resize(nVertices, nVertices);
	weight_.setZero();
    // iterating through each face to work with every single triangle
    for (int face = 0; face < nFaces; face++) {
        //loop over the three vertices within the face
        for (int v = 0; v < 3; v++) {
            //indices and coordinates of each vertex of each face
            firstVertex = faces_(face, VertexToEdge_map[v][0]);
            secondVertex = faces_(face, VertexToEdge_map[v][1]);
            thirdVertex = faces_(face, VertexToEdge_map[v][2]);

            A = vertices_.row(firstVertex);
            B = vertices_.row(secondVertex);
            C = vertices_.row(thirdVertex);
        }
        // get each angle in the face
        double alpha = angleBetweenVectors(A - B, A - C);
        double beta = angleBetweenVectors(B - A, B - C);
        double gamma = angleBetweenVectors(A - C, B - C);

        /*
			apply the weighting function to each vertex
			indexing is according to the paper
			and set s.t. the matrix is symmetric (weight_ = weight_.transpose
			*/
        // weights for the upper half
        weight_.coeffRef(secondVertex, thirdVertex) += 1 / (2.0 * std::tan(alpha));
        weight_.coeffRef(firstVertex, thirdVertex)  += 1 / (2.0 * std::tan(beta));
        weight_.coeffRef(firstVertex, secondVertex) += 1 / (2.0 * std::tan(gamma));
        // weights for the lower half
        weight_.coeffRef(thirdVertex, secondVertex) += 1 / (2.0 * std::tan(alpha));
        weight_.coeffRef(thirdVertex, firstVertex)  += 1 / (2.0 * std::tan(beta));
        weight_.coeffRef(secondVertex, firstVertex) += 1 / (2.0 * std::tan(gamma));
    }

    std::cout << "Computing weights ... DONE !" << std::endl;
}

double ArapCompute::angleBetweenVectors(const Eigen::Vector3d& a, const Eigen::Vector3d& b)
{
    double angle = 0.0;

    angle = std::atan2(a.cross(b).norm(), a.dot(b));

    return angle;
}

Eigen::VectorXi ArapCompute::computeNeighbourVertices(int vertexID){
    // returns a vectors Xi storing the IDs of all the neighbouring vertices to one vertex
    // in other words the cell related to one vertex
    // N(i) can be obtained by vertices.size()
    // TODO
    // find a more efficient solution
    int nVertices = vertices_.rows();
    int nFaces = faces_.size();

    Eigen::VectorXi neighbouringVertices;

    int firstVertex;
    int secondVertex;
    int thirdVertex;

    for (int face = 0; face < nFaces;)
    {
        for (int v = 0; v < 3; v++)
        {
            firstVertex  = faces_(face, VertexToEdge_map[v][0]);
            secondVertex = faces_(face, VertexToEdge_map[v][1]);
            thirdVertex  = faces_(face, VertexToEdge_map[v][2]);
        }      

        if (firstVertex == vertexID)
        {
            neighbouringVertices(neighbouringVertices.size() + 1) = secondVertex;
            neighbouringVertices(neighbouringVertices.size() + 1) = thirdVertex;
        }
        else if (secondVertex == vertexID)
        {
            neighbouringVertices(neighbouringVertices.size() + 1) = firstVertex;
            neighbouringVertices(neighbouringVertices.size() + 1) = thirdVertex;
        }
        else if (thirdVertex == vertexID)
        {
            neighbouringVertices(neighbouringVertices.size() + 1) = secondVertex;
            neighbouringVertices(neighbouringVertices.size() + 1) = thirdVertex;
        }
    }
    return neighbouringVertices;
}

void ArapCompute::computeLaplaceBeltramiOperator()
{
    std::cout << "Compute Laplace-Beltrami Operator ..." << std::endl;

    int nVertices = vertices_.rows();

	//compute the laplace-beltrami operator
	L_operator.resize(nVertices, nVertices);

	Eigen::VectorXi neighbours;
	double finalWeight;

	for (int i=0; i < nVertices; i++)
    {
        neighbours = computeNeighbourVertices(i);
        for (int j=0; i < neighbours.size(); j++){
        finalWeight = weight_.coeff(i,j);
        L_operator.coeffRef(i,i) += finalWeight;
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
    // computes the first guess for the vertex transformation $\pmb{p}_0'$ that is the basis for the
    // initial rotation
    std::cout << "Compute Initial guess ..." << std::endl;

    int n = vertices_.rows();

    Eigen::SparseMatrix<double> I(n, n);
    I.setIdentity();

    auto L_trans_L = L_operator.transpose() * L_operator;
    sparse_solver.compute(L_trans_L);

    auto L_Operator_inv = sparse_solver.solve(I);

    Eigen::VectorXd delta = L_operator * vertices_;
    // solution to $||Lp'-\delta||^2$: p' = (L.transpose * L).inverser()*L.transpose * delta
    vertices_ = L_Operator_inv * L_operator.transpose() * delta;

    std::cout << "Naive Laplacian Editing ... DONE !" << std::endl;
}

void ArapCompute::ComputeRotations()
{
    std::cout << "Compute rotations ..." << std::endl;

	//check equations (5) and (6) in the paper for computing rotations

	//total number of vertices
	int nVertices = vertices_.rows();

    Eigen::VectorXd neighbours;
    Eigen::Vector3d restEdge;
    Eigen::Vector3d deformedEdge;

	//a vector of matrices to hold the products in equation (5) in the paper.
	//S = per-edge weights * rest edge * deformed edge, summed over all vertices.
	std::vector<Eigen::MatrixXd> S_matrix;

	for (int v_1 = 0; v_1 < nVertices; v_1++)
	{
        neighbours = computeNeighbourVertices(v_1);
		for (int v_2 = 0; v_2 < neighbours.size(); v_2++)
		{
			restEdge     = vertices_.row(v_1) - vertices_.row(v_2);
			deformedEdge = updatedVertices_.row(v_1) - updatedVertices_.row(v_2);

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


    for (int i = 0; i < nVertices; i++)
    {
        Eigen::Vector3f sum(0.0f, 0.0f, 0.0f);
        for (int idx = 0;idx < neighbours.size(); idx++)
        {
            unsigned int j = neighbours(idx);
            sum += weight_.coeff(i,j) / 2.0f * (rotations[i] + rotations[j]) * (vertices_[i] - vertices_[j]);
        }

        RHS(i,0) = sum.x();
        RHS(i,1) = sum.y();
        RHS(i,2) = sum.z();
    }

    std::cout << "Compute right hand side ... DONE!" << std::endl;
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
