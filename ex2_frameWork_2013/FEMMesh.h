/*=====================================================================================*/
/*! \file		FEMMesh.h
	\author		peterkau
	\brief		Declaration of class FEMMesh
 */
/*=====================================================================================*/

#ifndef FEM_MESH_H
#define FEM_MESH_H

#include "FEMElementTri.h"
//#include "IFEMAssembler.h"

//! Nodal boundary condition (Dirichlet boundary condition)
class BoundaryCondition
{
public:
	BoundaryCondition(size_t nodeID, float val) {
		m_nodeID = nodeID;
		m_val = val;
	}

	//! Gets the global ID of the node at which the boundary condition applies
	size_t GetID() const { return m_nodeID; }

	//! Gets the value of the Dirichlet boundary condition
	float GetValue() const { return m_val; }

private:
	//! Node at which the boundary condition applies
	size_t m_nodeID;

	//! Value of the Dirichlet boundary condition
	float m_val;
};

//! A mesh consisting of nodes and elements
class FEMMesh {
public:
	//! Adds a node to the mesh and returns the ID of the node
	size_t AddNode(const Vector2 &position) {
		m_nodes.push_back(position);
		return m_nodes.size() - 1;
	}

	//! Adds an element to the mesh and returns the ID of the element
	size_t AddElement(FEMElementTri element) {
		m_elements.push_back(element);
		return m_elements.size() - 1;
	}

	//! Sets the right-hand side of the problem
	void SetRHS(const vector<float> &rhs) { m_rhs = rhs; }

	//! Returns the number of nodes (= degrees of freedom) in the mesh
	size_t GetNumNodes() const {
		return m_nodes.size();
	}

	//! Returns the number of elements in the mesh
	size_t GetNumElements() const {
		return m_elements.size();
	}

	//! Returns the element with ID \c elID
	const FEMElementTri &GetElement(size_t elID) const {
		return m_elements[elID];
	}

	// From IFEMAssembler:
	virtual const Vector2 &GetNodePosition(size_t nodeID) const {
		return m_nodes[nodeID];
	}

	// From IFEMAssembler:
	virtual void AddToStiffnessMatrix(size_t u, size_t v, float val) {
		m_matK(u, v) += val;
	}

	//! access Stiffness matrix
	SparseSymmetricDynamicRowMatrix& getMat() { return m_matK; }


	//! Solves the problem. NOTE: Call \c BuildStiffnessMatrix() before calling this
	bool Solve(vector<float> &solution, const vector<BoundaryCondition> &boundaryConds) const
	{
		if(m_matK.GetNumRows() != m_rhs.size()) {
			std::cerr<<"FEMMesh::Solve matrix and RHS dont match! "<<m_matK.GetNumRows()<<" vs "<<m_rhs.size() <<" \n";
			return false;	// Dimension mismatch! RHS not computed?
		}

		solution.resize(m_matK.GetNumRows());

		// Create copies of K and RHS:
		vector<float> rhs = m_rhs;
		SparseSymmetricDynamicRowMatrix matK = m_matK;

		// Apply Dirichlet boundary conditions by modifying K and RHS:
		for(size_t i=0; i<boundaryConds.size(); i++)
		{
			matK.FixSolution(rhs, boundaryConds[i].GetID(), boundaryConds[i].GetValue());
		}

		// Solve the problem:
		SparseLinSolverPCGT<float> solver;
		return solver.SolveLinearSystem(matK, solution, rhs, (float)1e-6, 1000);
	}

	void printStiffnessMatrix() const
	{
		int n = m_matK.GetNumRows();
		for(int i =0; i<n;i++)
		{
			printf("[");
			for(int j=0; j<n;j++)
			{
				double d = (i>j)?m_matK(i,j):m_matK(j,i);
				printf("%f ",d);
			}
			printf("]\n");
		}
	}

private:
	//! Position of all nodes in the mesh
	vector<Vector2> m_nodes;

	//! All elements in the mesh
	vector<FEMElementTri> m_elements;

	//! Stiffness matrix
	SparseSymmetricDynamicRowMatrix m_matK;

	//! Current right-hand side
	vector<float> m_rhs;
};

#endif
