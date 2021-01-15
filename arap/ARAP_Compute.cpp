#include "ARAP_Compute.h"

#include "igl/slice.h"

//Helps mapping from a vertex to the corresponding edge
const int VertexToEdge_map[3][2] = { {0, 1}, {1, 2}, {2, 0} };


ArapCompute::ArapCompute(const Eigen::MatrixXd& vertices, const Eigen::VectorXi& fixedVertices,
	const Eigen::MatrixXi& faces, const int maxIterations)
	:Compute(vertices, fixedVertices, faces, maxIterations) {}



void ArapCompute::Prefactorize()
{
	int nVertices = vertices_.rows();
	int nFree = freeVertices_Index.size();
	int nFaces = faces_.size();

	//compute weights
	weight_.resize(nVertices, nVertices);

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

			//using .coeff() gave an error but using coeffRef() didn't
			//in Eigen::SparseMatrix, coeffRef() returns a non-constant reference to the value
			//of the matrix at position i, j
			weight_.coeffRef(firstVertix, secondVertix) += cotan(v) / 2.0;

			//since the weight of the edge i-j (edge between vertix(i) and vertix(j) is the same
			//of the edge j-i, then we assign them the same weight in the weight_ matrix
			weight_.coeffRef(secondVertix, firstVertix) += cotan(v) / 2.0;
		}
	}

	//compute neighbors
	neighbors_.resize(nVertices, vertixNeighbors());
	
	for (int face = 0; face < nFaces; face++)
	{
		for (int v = 0; v < nVertices; v++)
		{
			int firstVertix = faces_(face, VertexToEdge_map[v][0]);
			int secondVertix = faces_(face, VertexToEdge_map[v][1]);

			neighbors_[firstVertix][secondVertix] = secondVertix;
			neighbors_[secondVertix][firstVertix] = firstVertix;
		}
	}


	//compute the laplace-beltrami operator
	L_operator.resize(nFree, nFree);
	igl::slice(weight_, freeVertices_Index, freeVertices_Index, L_operator);
	L_operator *= -1.0;

	//for reducing memory consumption
	//makeCompressed() turns the sparseMatrix L_operator into the Compressed row/col storage form.
	L_operator.makeCompressed();

	//Sparse factorization
	sparse_solver.compute(L_operator);

	//TODO: Check for success of factorization
}

//Idea: Laplacian smoothing for the mesh 
//will leave it to the end.
void ArapCompute::Preprocess(const Eigen::MatrixXd& fixedVertices)
{
	fixedVertices_ = fixedVertices;

	int nVertices = vertices_.rows();
	int nFaces = faces_.size();
	int nFree = freeVertices_Index.size();


}

void ArapCompute::SingleIteration()
{
	int nVertices = vertices_.rows();
	int nFixedVertices = fixedVertices_Index.size();
	int nFaces = faces_.size();

	//a vector of matrices to hold the products in equation (5) in the paper.
	//S = the sum over all vertices (per-edge weights * rest edge * deformed edge)
	std::vector<Eigen::MatrixXd> S_matrix;

	//Step 1: Estimate rotations using SVD
	ComputeRotations();


	//Step 2: Compute the matrix of rotated vertex gradients (the R.H.S in equation (9))
	// equation (9): L * P' = B, where B is an nx3 matrix
	// B = KR 
	// K is nx3n sparse matrix containing Wij * (Pi - Pj)
	// R is a matrix of rotation matrices
	//Reminder: n is the number of vertices in the mesh



}

Eigen::Vector3d ArapCompute::ComputeCotangent(int face_index) const
{
	//inititalize the cotangent vector
	//REMINDER: the cotangent is the per-edge weight in the paper (i.e. Wij)
	Eigen::Vector3d cotangent(0.0, 0.0, 0.0);

	return cotangent;
}

double ArapCompute::ComputeEnergy() const
{

}

std::Vector<double> ComputeFaceArea(const Eigen::MatrixXi& faces) const
{

}




