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

void AdvanceTimeStep1(double k, double m, double d, double L, double dt, int method, double p1, double v1, double& p2, double& v2)
{
	double x, v;

	switch(method){
	 // euler
	 case 1:
		x = p2 + dt*v2;
		v = v2 + dt*((-k*(L+p2) - d*v2)/m - g);

		p2 = x;
		v2 = v;
		break;
	 // semi-implicit euler / symplectic euler
	 case 2:
		v = v2 + dt*((-k*(L+p2) - d*v2)/m - g);
		x = p2 + dt*v;

		p2 = x;
		v2 = v;
		break;
	 // midpoint
	 case 3:
		 double vMid, xMid;
		 xMid = p2 + dt/2*v2;
		 vMid = v2 + dt/2*((-k*(L+p2) - d*v2)/m - g);

		 x = p2 + dt*vMid;
		 v = vMid + dt*((-k*(L+xMid) - d*vMid)/m - g);

		 p2 = x;
		 v2 = v;
	 // semi-implicit backward euler (to be implemented)

	 // analytic
	 case 5:
		 double c1, alpha, beta;
		 time += dt;
		 c1 = m*g/k;
		 alpha = -d/(2*m);
		 beta = sqrt(4*k*m - d*d)/(2*m);

		 x = c1*exp(alpha*time)*cos(beta*time) - alpha/beta*c1*exp(alpha*time)*sin(beta*time) - L - (m*g)/k;
		 v = -c1*(alpha*alpha + beta*beta)/beta*sin(beta*time);

		 p2 = x;
		 v2 = v;

	}

}


// Exercise 3
// falling triangle
void AdvanceTimeStep3(double k, double m, double d, double L, double dt,
                      Vec2& p1, Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3)
{
   p1 += Vec2(1,1);
}
