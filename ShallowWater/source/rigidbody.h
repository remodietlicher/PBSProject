#ifndef RIGIDBODY_H
#define RIGIDBODY_H

#include "utilities.h"

struct RigidBody{
	// constants
	float mass;
	Matrix3x3f Ibody, Ibodyinv;

	// state
	Vector3f x;
	Matrix3x3f R;
	Vector3f P, L;

	// auxiliary
	Matrix3x3f Iinv;
	Vector3f v, omega;

	// computed quantities
	Vector3f force, torque;
	RigidBody(Vector3f com, float m) :
		mass(m),
		x(com),
		R(1, 0, 0, 0, 1, 0, 0, 0, 1), // Identity
		P(0, 0, 0),
		L(0, 0, 0),
		v(0.0f, 0.0f, 0.0f)
	{}
};

struct Box : RigidBody{
	Vector3f x0, y0, z0; // axes in x, y, z direction
	Box(Vector3f com, float m, Vector3f a, Vector3f b, Vector3f c) :
		RigidBody(com, m),
		x0(a),
		y0(b),
		z0(c)
	{
		// inertia tensor of a box
		Ibody = mass/12.0f*Matrix3x3f(y0.squaredLength()+z0.squaredLength(), 0.0f, 0.0f, 0.0f, x0.squaredLength()+z0.squaredLength(), 0.0f, 0.0f, 0.0f, x0.squaredLength()+y0.squaredLength());
		Ibodyinv = Matrix3x3f(1.0f/Ibody(0, 0), 0, 0, 0, 1.0f/Ibody(1, 1), 0, 0, 0, 1.0f/Ibody(2, 2));
//		Ibodyinv = Ibody.inverse();
		Iinv = R*Ibodyinv*R.transposed();
		omega = Iinv*L;
	}
};

#endif