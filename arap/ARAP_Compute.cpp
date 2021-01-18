#include "ARAP_Compute.h"

#include "igl/slice.h"
#include "igl/polar_svd.h"
#include "igl/polar_svd3x3.h"


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
    //compute weights
    int nVertices = vertices_.rows();
    int nFaces = faces_.size();

	weight_.resize(nVertices, nVertices);
	weight_.setZero();
    // iterating through each face to work with every single triangle
    for (int face = 0; face < nFaces; face++) {
        //loop over the three vertices within the face
        for (int v = 0; v < 3; v++) {
            //indices and coordinates of each vertex of each face
            int firstVertex = faces_(face, VertexToEdge_map[v][0]);
            int secondVertex = faces_(face, VertexToEdge_map[v][1]);
            int thirdVertex = faces_(face, VertexToEdge_map[v][2]);

            Eigen::Vector3d A = vertices_.row(firstVertex);
            Eigen::Vector3d B = vertices_.row(secondVertex);
            Eigen::Vector3d C = vertices_.row(thirdVertex);
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

    for (int face = 0; face < nFaces;)
    {
        int firstVertex  = faces_(face, VertexToEdge_map[v][0]);
		int secondVertex = faces_(face, VertexToEdge_map[v][1]);
        int thirdVertex  = faces_(face, VertexToEdge_map[v][2]);

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

    int nVertices = vertices_.rows();

	//compute the laplace-beltrami operator
	L_operator.resize(nVertices, nVertices);

	Eigen::VectorXi neighbours;
	double finalWeight;

	for (int i=0; i < nVertices; i++){
        neighbours = computeNeighbourVertices(i);
        for (int j=0; i < neighbours.size(); j++){
        finalWeight = weight_.coeff(i,j);
        L_operator(i,i) += finalWeight;
        }
	}

	//for reducing memory consumption
	//makeCompressed() turns the sparseMatrix L_operator into the Compressed row/col storage form.
	L_operator.makeCompressed();
}

void ArapCompute::NaiveLaplaceEditing() {
    // computes the first guess for the vertex transformation $\pmb{p}_0'$ that is the basis for the
    // initial rotation

    Eigen::VectorXd delta = L_operator * vertices_;
    // solution to $||Lp'-\delta||^2$: p' = (L.transpose * L).inverser()*L.transpose * delta
    vertices_ = (L_operator.transpose() * L_operator).inverse() * L_operator.transpose() * delta;

}

void ArapCompute::ComputeRotations()
{
    VERBOSE("Compute rotations ...");
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
			restEdge = vertices_.row(v_1) - vertices_.row(v_2);
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

		igl::polar_svd3x3(S_matrix[v], rotation);
		rotations[v] = rotation.transpose();
	}
    VERBOSE("Compute rotations ... DONE!");
}

void ArapCompute::ComputeRightHandSide()
{
    VERBOSE("Compute right hand side ...");



    Eigen::VectorXi neighbours;


    for (int i = 0; i < (int)m_nVertices; ++i)
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

    VERBOSE("Compute right hand side ... DONE!");
}


void ArapCompute::alternatingOptimization()
{
    for (int iter = 0; iter < maxIterations_; iter++)
    {

    }
    // update the vertices


}
