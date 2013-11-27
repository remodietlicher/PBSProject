#ifndef RBS_H
#define RBS_H

#include "rigidbody.h"

class RigidBodySolver{
private:
	Box body;
public:
	RigidBodySolver(Box body);
	void advanceTimestep(float dt, Vector3f* F, Vector3f* r, int N);
	Box* getBody();

private:
	Matrix3x3f star(Vector3f v);
	void sumExternalForces(Vector3f* F, Vector3f* r, int N);
	void performEulerIntegration(float dt);
	void setBodyBox();
};

#endif