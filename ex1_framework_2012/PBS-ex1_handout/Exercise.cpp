#include <iostream>
using namespace std;
#include "Vec2.h"
#include "LinSys.h"
#include <math.h>

// gravitational acceleration (9.81)
static const double g = 9.81f;


// Exercise 1
// hanging mass point
// param k : spring constant
// param m : mass of p2
// param d : damping term
// param dt: timestep
// param k : spring constant

static double time = 0;

Vec2 calculateForce2D(Vec2 x0, Vec2 xi, Vec2 v, bool isOnGround, double k, double m, double d, double L);
double calculateForce1D(double k, double m, double d, double L, double p2, double v2);

void AdvanceTimeStep1(double k, double m, double d, double L, double dt, int method, double p1, double v1, double& p2, double& v2)
{
	double x, v;

	switch(method){
	 // euler
	 case 1:
		x = p2 + dt*v2;
		v = v2 + dt*calculateForce1D(k,m,d,L,p2,v2);
		break;
	 // semi-implicit euler / symplectic euler
	 case 2:
		v = v2 + dt*calculateForce1D(k,m,d,L,p2,v2);
		x = p2 + dt*v;
		break;
	 // midpoint
	 case 3:
		 double vMid, xMid;
		 xMid = p2 + dt/2*v2;
		 vMid = v2 + dt/2*calculateForce1D(k,m,d,L,p2,v2);

		 x = p2 + dt*vMid;
		 v = vMid + dt*calculateForce1D(k,m,d,L,xMid,vMid);
		 break;
	 // semi-implicit backward euler (TODO: to be implemented)


	 // analytic
	 case 5:
		 double c1, alpha, beta;
		 time += dt;
		 c1 = m*g/k;
		 alpha = -d/(2*m);
		 beta = sqrt(4*k*m - d*d)/(2*m);

		 x = c1*exp(alpha*time)*cos(beta*time) - alpha/beta*c1*exp(alpha*time)*sin(beta*time) - L - (m*g)/k;
		 v = -c1*(alpha*alpha + beta*beta)/beta*sin(beta*time);
		 break;
	}

	p2 = x;
	v2 = v;
}


// Exercise 3
// falling triangle
void AdvanceTimeStep3(double k, double m, double d, double L, double dt,
                      Vec2& p1, Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3)
{
	Vec2* points[3] = {&p1, &p2, &p3};
	Vec2* velocities[3] = {&v1, &v2, &v3};
	
	// symplectic euler
	Vec2 v, x, v_1, v_2;

	// looping through all masspoints
	for(int i=0; i<3; i++){
		Vec2* pointA = points[i];
		Vec2* pointB = points[(i+1)%3];
		Vec2* pointC = points[(i+2)%3];
		Vec2* vA = velocities[i];

		// check if the mass under consideration (mass at pointA)
		bool isOnGround;
		if((*pointA).y <= -1){
			isOnGround = true;
		} else isOnGround = false;

		// velocity step (t --> t+dt) for interaction between pointA and pointB
		v_1 = dt*calculateForce2D(*pointA, *pointB, *vA, isOnGround, k, m, d, L);
		// velocity step for interaction between pointA and pointC
		v_2 = dt*calculateForce2D(*pointA, *pointC, *vA, isOnGround, k, m, d, L);
		// total velocity step
		v = *vA + v_1 + v_2;
		// coordinate step (t --> t + dt) 
		x = *pointA + dt*v;

		//if(x.y < -1){
		//	x.y = (*pointA).y;
		//}

		*pointA = x;
		*vA = v;
	}
}

// calculate the force due to a spring of length L between x and the origin
double calculateForce1D(double k, double m, double d, double L, double x, double v){
	return (-k*(L+x) - d*v)/m - g;
}

// calculate the force due to a spring of length L between x0 and xi
Vec2 calculateForce2D(Vec2 x0, Vec2 xi, Vec2 v, bool isOnGround, double k, double m, double d, double L){

	Vec2 Fspring, Fdamp, Fgrav, Fground, Ftot;
	double kground = 200.0;

	// Force due to Springs (discarded "-" in front of "k" to get a force in the right direction e.g.
	// if |xi - x0| > L (streched spring) we want that x0 gets pulled towards xi)
	Fspring = k*((xi-x0).length()-L) * 1.0/(xi-x0).length()*(xi-x0);
	// Damping Force
	Fdamp = -d*v;
	// Gravitational force
	Fgrav = g*Vec2(0, -1.0);

	Ftot = Fspring + Fdamp + Fgrav;

	// handle ground-interaction
	if(isOnGround) {
		Fground = -kground*(-1.0-x0.y)*Vec2(0, 1);
		Ftot -= Fground;
	}

	Ftot = 1.0/m * Ftot;

	return Ftot;
}
