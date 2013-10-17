/*=====================================================================================*/
/*! \file		FEMElementTri.h
	\author		peterkau
	\brief		Declaration of class FEMElementTri
 */
/*=====================================================================================*/

#ifndef FEM_ELEMENT_TRI_H
#define FEM_ELEMENT_TRI_H

class FEMMesh;

class FEMElementTri
{
public:
	FEMElementTri(size_t node0, size_t node1, size_t node2) {
		m_nodes[0] = node0;
		m_nodes[1] = node1;
		m_nodes[2] = node2;
	}

	//! Adds the contributions of this element to the global stiffness matrix
	void Assemble(FEMMesh *pMesh) const;

	//! Returns the number of nodes in the element
	size_t GetNumElementNodes() const { return 3; }

	//! Returns the global ID of the element node \c elsize_t
	size_t GetGlobalNodeForElementNode(size_t elsize_t) const { return m_nodes[elsize_t]; }

//private:

	/*! Computes for a basis function its x/y derivatives geometrically*/
	void computeSingleBasisDerivGlobalGeom(size_t nodeId, Vector2 &basisDerivGlobal, const FEMMesh *pMesh) const;

	/*! Computes for a basis function its x/y derivatives geometrically*/
	void computeSingleBasisDerivGlobalLES(size_t nodeId, Vector2 &basisDerivGlobal, const FEMMesh *pMesh) const;

private:
	//! The global IDs of the three nodes of the triangle
	size_t m_nodes[3];
};

#endif
