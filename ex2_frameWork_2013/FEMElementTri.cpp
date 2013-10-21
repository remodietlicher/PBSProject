/*=====================================================================================*/
/*! \file		FEMElementTri.cpp
	\author		peterkau
	\brief		Implementation of class FEMElementTri
 */
/*=====================================================================================*/

#include "SimpleFEMDefs.h"
#include "FEMElementTri.h"
#include "FEMMesh.h"
#include <iostream>
#include <math.h>

// choose between geometric and usual gradient computation
#define computeSingleBasisDerivGlobalVar computeSingleBasisDerivGlobalLES

float FEMElementTri::getAreaInMesh(const FEMMesh mesh) const {
	Vector2 x0 = mesh.GetNodePosition(GetGlobalNodeForElementNode(0));
	Vector2 x1 = mesh.GetNodePosition(GetGlobalNodeForElementNode(1));
	Vector2 x2 = mesh.GetNodePosition(GetGlobalNodeForElementNode(2));
	Vector2 x12 = x2-x1;
	Vector2 x10 = x0-x1;
	float area = 0.5*abs((x12[0]*x10[1]-x12[1]*x10[0]));
	
	return area;
}

Vector3 FEMElementTri::calculateCoefficientsInMesh(const FEMMesh mesh, size_t nodeId) const {
	Vector2 x0 = mesh.GetNodePosition(GetGlobalNodeForElementNode(0));
	Vector2 x1 = mesh.GetNodePosition(GetGlobalNodeForElementNode(1));
	Vector2 x2 = mesh.GetNodePosition(GetGlobalNodeForElementNode(2));
	Matrix3x3 A(x0[0], x0[1], 1, x1[0], x1[1], 1, x2[0], x2[1], 1);

	// note: nodeId is given in the element-system, therefore 0<nodeId<2
	// build vector b = (delta_1i, delta_2i, delta_3i)
	assert(nodeId<3);
	Vector3 b = Vector3(0, 0, 0);
	b[nodeId] = 1;
	return A.inverse()*b;
}

float FEMElementTri::eval_N(const FEMMesh mesh, size_t nodeId, Vector2 x) const {
	Vector3 coeffs = calculateCoefficientsInMesh(mesh, nodeId);
	return x[0]*coeffs[0] + x[1]*coeffs[1] + coeffs[2];
}

// TASK 3
void FEMElementTri::Assemble(FEMMesh *pMesh) const
{
	//float area = getAreaInMesh(*pMesh); --> wieso ist das so langsam?

//	/*
	Vector2 x0 = pMesh->GetNodePosition(GetGlobalNodeForElementNode(0));
	Vector2 x1 = pMesh->GetNodePosition(GetGlobalNodeForElementNode(1));
	Vector2 x2 = pMesh->GetNodePosition(GetGlobalNodeForElementNode(2));
	Vector2 x12 = x2-x1;
	Vector2 x10 = x0-x1;
	float area = 0.5*abs((x12[0]*x10[1]-x12[1]*x10[0]));
//	*/

	Vector2 basisDerivGlobal_i, basisDerivGlobal_j;
	int global_i, global_j;
	float value;
	for(int i=0; i<3; i++){
		computeSingleBasisDerivGlobalVar(i, basisDerivGlobal_i, pMesh);
		for(int j=0; j<=i; j++){
			computeSingleBasisDerivGlobalVar(j, basisDerivGlobal_j, pMesh);
			value = area*(basisDerivGlobal_i|basisDerivGlobal_j);
			// can only write to the lower part of K_ij, so we must assert that Global_i (row) >= Global_j (column)
			global_i = GetGlobalNodeForElementNode(i);
			global_j = GetGlobalNodeForElementNode(j);
			if(global_i < global_j){
				int tmp = global_i;
				global_i = global_j;
				global_j = tmp;
			}
			pMesh->AddToStiffnessMatrix(global_i, global_j, value);
		}
	}
}

// TASK 2
void FEMElementTri::computeSingleBasisDerivGlobalGeom(size_t nodeId, Vector2 &basisDerivGlobal, const FEMMesh *pMesh) const
{
	Vector2 x0 = pMesh->GetNodePosition(GetGlobalNodeForElementNode(nodeId));
	Vector2 x1 = pMesh->GetNodePosition(GetGlobalNodeForElementNode((nodeId+1)%3));
	Vector2 x2 = pMesh->GetNodePosition(GetGlobalNodeForElementNode((nodeId+2)%3));
	Vector2 x12 = x2-x1;
	Vector2 x10 = x0-x1;

	// calculate normal vector to x12
	Vector2 nx12 = Vector2(-x12[1], x12[0]);
	assert(nx12[0]*x12[0] + nx12[1]*x12[1] == 0);
	float height = abs((x12[0]*x10[1]-x12[1]*x10[0]))/x12.length(); // height = 2*area/base
	basisDerivGlobal = nx12/nx12.length()*(1/height); // vector normal to x12 of unit (=inverse of triangles height) length
}

// TASK 1
void FEMElementTri::computeSingleBasisDerivGlobalLES(size_t nodeId, Vector2 &basisDerivGlobal, const FEMMesh *pMesh) const
{
	//Vector3 coeffs = calculateCoefficientsInMesh(*pMesh, nodeId); --> wieso ist das so langsam?

//	/*
	Vector2 x0 = pMesh->GetNodePosition(GetGlobalNodeForElementNode(0));
	Vector2 x1 = pMesh->GetNodePosition(GetGlobalNodeForElementNode(1));
	Vector2 x2 = pMesh->GetNodePosition(GetGlobalNodeForElementNode(2));
	Matrix3x3 A(x0[0], x0[1], 1, x1[0], x1[1], 1, x2[0], x2[1], 1);

	// note: nodeId is given in the element-system, therefore 0<nodeId<2
	// build vector b = (delta_1i, delta_2i, delta_3i)
	assert(nodeId<3);
	Vector3 b = Vector3(0, 0, 0);
	b[nodeId] = 1;
	Vector3 coeffs = A.inverse()*b;
//	*/

	basisDerivGlobal = Vector2(coeffs[0], coeffs[1]);
}
