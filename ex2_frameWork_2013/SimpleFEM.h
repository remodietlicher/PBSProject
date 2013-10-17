/*=====================================================================================*/
/*! \file		SimpleFEM.h
	\author		peterkau
	\brief		Declaration of class SimpleFEM
 */
/*=====================================================================================*/

#ifndef SIMPLE_FEM_H
#define SIMPLE_FEM_H

#include "FEMMesh.h"


//! Class containing various helper functions to test the FEMMesh class
class SimpleFEM
{
public:
	/*! Creates a uniform FEM triangle mesh with \c nodesX * \c nodesY nodes for the
		domain [0,0] - [1,1]. */
	static FEMMesh CreateUniformGridMesh(size_t nodesX, size_t nodesY);

	/*! Computes the boundary conditions for a known solution on the domain [0,0] - [1,1]. */
	static void ComputeBoundaryConditions(const FEMMesh &mesh, vector<BoundaryCondition> &boundaryConditions);
	/*! Computes the right-hand side for a known solution on the domain [0,0] - [1,1]. */
	static void ComputeRHS(const FEMMesh &mesh, vector<float> &rhs);
	/*! Computes the error of a numerical solution and its energy norm */
	static void computeError(FEMMesh &mesh,  const vector<float> &sol_num, vector<float> &verror, float& err_nrm);

	//! Returns \c true iff the point \c pos is on the boundary of the domain [0,0] - [1,1]
	static bool isOnBoundary(const Vector2 &pos);
};

#endif
