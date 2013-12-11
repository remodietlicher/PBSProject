#ifndef RBS_H
#define RBS_H

#include "rigidbody.h"

class RigidBodySolver{
private:
	Box* body;
public:
	RigidBodySolver(Box *body);
	void advanceTimestep(float dt, std::vector<Vector3f> F, std::vector<Vector3f> r);
	Box* getBody();

private:
	Matrix3x3f star(Vector3f v);
	void sumExternalForces(std::vector<Vector3f> F, std::vector<Vector3f> r);
	void performEulerIntegration(float dt);
	void setBodyBox();
};

#endif