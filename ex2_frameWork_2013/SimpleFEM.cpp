#include "SimpleFEMDefs.h"
#include "SimpleFEM.h"
#include "MeshViewer.h"

// size of grid
#ifndef GRIDSIZE // you can set GRIDSIZE with compiler flag -DGRIDSIZE=...
	static const size_t GRIDSIZE = 16;
#endif
// use a graded mesh, or a regular mesh
static const bool gradedMesh = false;
// laplace or poisson problem?
static const bool laplaceProblem = false;
// plot solution or error?
static bool vizSolution = false;
// display debug information?
static const bool debugOut = false;


float eval_u(float x, float y)
{
	if(laplaceProblem)
		//return x*x - y*y;
		return exp(x)*sin(y);
	else
		return 	3.f*x*x + 2.f*y*y*y*x;
}

float eval_f(float x, float y)
{
	if(laplaceProblem)
		return 0;
	else
		return 	-6.f - 12.f*y*x;
}


FEMMesh SimpleFEM::CreateUniformGridMesh(size_t nodesX, size_t nodesY)
{
	assert(nodesX >= 2);
	assert(nodesY >= 2);

	FEMMesh mesh;

	// setup positions of nodes
	int nodecnt =0;
	for(size_t y=0; y<nodesY; y++)
	{
		for(size_t x=0; x<nodesX; x++)
		{
			Vector2 pos = Vector2((float)x / (float)(nodesX-1), (float)y / (float)(nodesY-1));

			//shift mesh positions for graded mesh
			if(gradedMesh)
			{
				pos[0] *= pos[0]; pos[1] *= pos[1];
				//pos[0] = sqrt(pos[0]); pos[1] = sqrt(pos[1]);
			}
			mesh.AddNode(pos);
			nodecnt++;
		}
	}
	std::cout<<"Added "<<nodecnt<<" nodes to mesh.\n";

	// create elements
	int cnt=0;
	for(size_t y=0; y<nodesY-1; y++)
	{
		for(size_t x=0; x<nodesX-1; x++)
		{
			// bottom-left:
			size_t node00 = y*nodesX + x;

			// bottom-right:
			size_t node10 = node00 + 1;

			// top-left:
			size_t node01 = node00 + nodesX;

			// top-right:
			size_t node11 = node00 + nodesX + 1;

			// add two element for this quad
			mesh.AddElement( FEMElementTri(node00, node10, node11) );
			mesh.AddElement( FEMElementTri(node00, node11, node01) );
			cnt+=2;
		}
	}
	std::cout<<"Added "<<cnt<<" elements to mesh.\n";

	return mesh;
}

void SimpleFEM::ComputeBoundaryConditions(const FEMMesh &mesh, vector<BoundaryCondition> &boundaryConditions)
{
	boundaryConditions.clear();

	for(size_t nodeID=0; nodeID<mesh.GetNumNodes(); nodeID++)
	{
		const Vector2 &pos = mesh.GetNodePosition(nodeID);

		if(isOnBoundary(pos)) {
			float x = pos[0];
			float y = pos[1];

			// compute reference solution on boundary
			//float val = 3.*x*x + 2.*y*y*y*x;
			float val = eval_u(x,y);

			// this fixes the solution for node "nodeID" to "val" when
			// solving the system later on
			boundaryConditions.push_back( BoundaryCondition(nodeID, val) );
		}
	}
}


// perform a simple boundary check
// is either of the components 0 or 1?
bool SimpleFEM::isOnBoundary(const Vector2 &pos)
{
	return pos[0] <= 0. || pos[0] >= 1. || pos[1] <= 0. || pos[1] >= 1.;
}


void SimpleFEM::ComputeRHS(const FEMMesh &mesh,  vector<float> &rhs)
{
	for(size_t ie=0; ie<mesh.GetNumElements(); ie++) {
		const FEMElementTri& elem = mesh.GetElement(ie);

		Vector2 x0 = mesh.GetNodePosition(elem.GetGlobalNodeForElementNode(0));
		Vector2 x1 = mesh.GetNodePosition(elem.GetGlobalNodeForElementNode(1));
		Vector2 x2 = mesh.GetNodePosition(elem.GetGlobalNodeForElementNode(2));
		Vector2 x12 = x2-x1;
		Vector2 x10 = x0-x1;
		float area = 0.5*fabs((x12[0]*x10[1]-x12[1]*x10[0]));

		Vector2 q = Vector2(0, 0);
		for(int i=0; i<3; i++){
			q += mesh.GetNodePosition(elem.GetGlobalNodeForElementNode(i));
		}
		q *= 1.0/3.0;

		Matrix3x3 A;
		for (int j = 0; j < 2; ++j)
		{
			A(0,j) = x0[j];
			A(1,j) = x1[j];
			A(2,j) = x2[j];
		}
		for (int j = 0; j < 3; ++j)
			A(j,2) = 1.0;

		float value;
		for(int i=0; i<3; i++){

			Vector3 b = Vector3(0, 0, 0);
			b[i] = 1;
			Vector3 coeffs = A.inverse()*b;

			value = area*eval_f(q[0], q[1])*(coeffs[0]*q[0]+coeffs[1]*q[1]+coeffs[2]);
			rhs.at(elem.GetGlobalNodeForElementNode(i)) += value;
		}
	}
}

void SimpleFEM::computeError(FEMMesh &mesh,  const vector<float> &sol_num, vector<float> &verror, float& err_nrm )
{
    //Task 5 starts here

	int N = sol_num.size(); // length of vectors

	// "v_err = |sol_n - sol_a|"
	for (int i = 0; i < N; ++i)
	{
		Vector2 x = mesh.GetNodePosition(i);
		verror[i] = fabs( eval_u(x[0], x[1]) - sol_num[i] );
	}

	// "err = wurzel(v_err*K*v_err)"
	std::vector<float> K_verror(N);
	mesh.getMat().MultVector(verror, K_verror);
	float err_sum = 0.0;
	for (int i = 0; i < N; ++i)
		err_sum += verror[i] * K_verror[i];
	err_nrm = sqrt(err_sum);

    //Task 5 ends here
}

int main(int argc, char *argv[])
{
	// Create a uniform mesh:
	FEMMesh mesh = SimpleFEM::CreateUniformGridMesh(GRIDSIZE, GRIDSIZE);
	int nNodes = mesh.GetNumNodes();
	// Build its stiffness matrix:
	// loop over all elements, and compute their contributions
	// for the equations of their respective nodes
	mesh.getMat().ClearResize( mesh.GetNumNodes() );
	for(size_t i=0; i<mesh.GetNumElements(); i++) {
		if(debugOut) cout<<"Assembling "<<i<<"\n";
		mesh.GetElement(i).Assemble( &mesh );
	}

	// Compute boundary conditions and right-hand side:
	vector<BoundaryCondition> boundaryConditions;

	SimpleFEM::ComputeBoundaryConditions(mesh, boundaryConditions);

	// Apply right-hand side:
	vector<float> rhs(nNodes);
	SimpleFEM::ComputeRHS(mesh, rhs);
	mesh.SetRHS(rhs);

	// Solve the problem, this calls a preconditioned CG solver
	// for the sparse matrix with right hand side rhs
	// all nodes stored in "boundaryConditions" are fixed to certain values
	vector<float> solution;
	bool isSolved = mesh.Solve(solution, boundaryConditions);
	assert(isSolved);

	// print matrix for boundary nodes
	if(debugOut)
		for(size_t i=0; i<mesh.GetNumNodes(); i++) {
			const Vector2 & pi = mesh.GetNodePosition(i);
			if(SimpleFEM::isOnBoundary(pi))
				continue;
			for(size_t j=0; j<mesh.GetNumNodes(); j++) {
				const Vector2 & pj = mesh.GetNodePosition(j);
				if(SimpleFEM::isOnBoundary(pj))
					continue;
				if(j>i) cout<<mesh.getMat()(j,i)<<"\t";
				else    cout<<mesh.getMat()(i,j)<<"\t";
			}
			std::cout<<std::endl;
		//cout<<" rhs="<<rhs[i]<<" \n";
	} // debug output

	float err_nrm = 0;
	std::vector<float> verr(nNodes);
	SimpleFEM::computeError(mesh,solution,verr,err_nrm);
	printf("Error norm is %f\n",err_nrm);
	// Visualize the solution:
	// draw the triangles with colors according to solution
	// blue means zero, red means maxValue.
	// the default problem goes from 0-5 , for other problems,
	// adjust the maxValue parameter below (values <0, or >maxValue
	// are clamped for the display)
	MeshViewer viewer(argc, argv);

	if(vizSolution)
		viewer.VisualizeSolution(mesh, solution);
	else
		viewer.VisualizeError(mesh, verr);


	return 0;
}


