#include "ARAP_Compute.h"
#include <iostream>

#include "igl/slice.h"
#include "igl/polar_svd.h"
#include "igl/polar_svd3x3.h"

using namespace std;


ArapCompute::ArapCompute(Eigen::MatrixXd& vertices,
                         Eigen::VectorXi& fixedVertices,
	                     const Eigen::MatrixXi& faces, 
                         int maxIterations)
	: vertices_(vertices),
      faces_(faces),
      fixedVertices_Index(fixedVertices),
      maxIterations_(maxIterations){}

void ArapCompute::ComputeWeights()
{
    std::cout << "Computing weights ..." << std::endl;

    //compute weights
    int nVertices = vertices_.rows();
    int nFaces = faces_.size();

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
    int nFree = freeVertices_Index.size();

    //Initialize the right hand side matrix of equations (8) and (9).
    RHS = Eigen::MatrixXd::Zero(nFree, 3);

    for (int i = 0; i < nFree; i++)
    {
        Eigen::Vector3d sum(0.0, 0.0, 0.0);

        for (int idx = 0; idx < neighbours.size(); idx++)
        {
            int j = neighbours[idx];
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
    for (int v = 0; v < nVertices; v++)
    {
        Eigen::VectorXi neighbors = computeNeighbourVertices(v);

        //Loop over all neighbors
        for (int n = 0; n < neighbors.size(); n++) 
        {
            int m = neighbors[n];

            double weight = weight_.coeff(n, m);
            double edge_energy = 0.0;

            //See equation (7) in the paper
            Eigen::Vector3d updatedVertices_vector = (cacheVertices_.row(n) - cacheVertices_.row(m)).transpose();
            Eigen::Vector3d oldVertices_vector = (vertices_.row(n) - vertices_.row(m)).transpose();

            total_energy += weight * (updatedVertices_vector - rotations[n] * oldVertices_vector).squaredNorm();
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