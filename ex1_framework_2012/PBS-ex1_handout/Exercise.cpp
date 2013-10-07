#include <iostream>
#include <fstream>
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
// param L : initial spring lenght
// param dt: timestep
// param k : spring constant


double calculateForce1D(double k, double m, double d, double L, double x, double v);
Vec2 calculateForce2D(Vec2 x0, Vec2 xi, Vec2 v, bool isOnGround, double k, double m, double d, double L);

class Integrator{

private:
	double m_time;//total time

	double m_k;//spring constant
	double m_L;//spring length
	double m_m;//mass
	double m_d;//damping
	static const double m_g = g;//gravitational acceleration


	double m_c1, m_alpha, m_beta;//integration constant for the analytic solution

private:
		static Integrator* s_instance;
		int cnt = 0;
		ofstream compare_pos_file;
		ofstream compare_vel_file;
	
public:

	// not thread safe but we are single threaded
	static Integrator* getinstance(double k, double m, double d, double L){

		if (s_instance == nullptr) {
			s_instance = new Integrator(k,m,d,L);
		}

		return s_instance;
	}

	Integrator(double k, double m, double d, double L){
		
		m_time = 0;
		m_m = m;
		m_k = k;
		m_L = L;
		m_d = d;

		m_c1 = m*g / k;
		m_alpha = -d / (2 * m);
		m_beta = sqrt(4 * k*m - d*d) / (2 * m);

		compare_pos_file.open("integrated_positions.csv");
		compare_vel_file.open("integrated_velocities.csv");
	}

	void timestep(double dt) {
		m_time += dt;
	}

	void euler(double dt, double& p2, double& v2){
		
		double x, v;

		x = p2 + dt*v2;
		v = v2 + dt*calculateForce1D(m_k, m_m, m_d, m_L, p2, v2);

		p2 = x;
		v2 = v;
	}

	void sympletic_euler(double dt, double& p2, double& v2){

		double x, v;

		v = v2 + dt*calculateForce1D(m_k, m_m, m_d, m_L, p2, v2);
		x = p2 + dt*v;

		p2 = x;
		v2 = v;
	}

	void mid_point(double dt, double& p2, double& v2){

		double x, v;
		double vMid, xMid;

		xMid = p2 + dt / 2 * v2;
		vMid = v2 + dt / 2 * calculateForce1D(m_k, m_m, m_d, m_L, p2, v2);

		x = p2 + dt*vMid;
		v = vMid + dt*calculateForce1D(m_k, m_m, m_d, m_L, xMid, vMid);

		p2 = x;
		v2 = v;

		p2 = x;
		v2 = v;
	}

	void analytic(double dt, double& p2, double& v2){

		double x, v;

		x = m_c1*exp(m_alpha*m_time)*cos(m_beta*m_time) - m_alpha / m_beta*m_c1*exp(m_alpha*m_time)*sin(m_beta*m_time) - m_L - (m_m*m_g) / m_k;
		v = -m_c1*(m_alpha*m_alpha + m_beta*m_beta) / m_beta*sin(m_beta*m_time);

		p2 = x;
		v2 = v;

	}


	void trace_all_methods(double dt, double p2, double v2) {

		if (++cnt > 1000) {
			return;
		}

		void(Integrator::* analytic) (double, double&, double&) = &Integrator::analytic;
		void(Integrator::* euler) (double, double&, double&) = &Integrator::euler;
		void(Integrator::* sympletic_euler) (double, double&, double&) = &Integrator::sympletic_euler;
		void(Integrator::* mid_point) (double, double&, double&) = &Integrator::mid_point;


		trace_method(analytic, dt, p2, v2);
		trace_method(euler, dt, p2, v2);
		trace_method(sympletic_euler, dt, p2, v2);
		trace_method(mid_point, dt, p2, v2);

		compare_pos_file << endl;
		compare_vel_file << endl;

		compare_pos_file.flush();
		compare_vel_file.flush();
	}

	void trace_method(void(Integrator::* method_ptr) (double, double&, double&), double dt, double p2, double v2) {

		(this->*method_ptr)(dt, p2, v2);

		compare_pos_file << p2 << ";";
		compare_vel_file << v2 << ";";
	}
};


//initialize global ptr
Integrator* Integrator::s_instance = nullptr;


void AdvanceTimeStep1(double k, double m, double d, double L, double dt, int method, double p1, double v1, double& p2, double& v2)
{
	// copy input references
	double in_p2 = p2, in_v2 = v2;

	// get singelton
	Integrator* integrator = Integrator::getinstance(k,m,d,L);

	// update total time
	integrator->timestep(dt);

	//trace all integrations
	integrator->trace_all_methods(dt, in_p2, in_v2);
	

	double x, v;

	switch(method){
	 // euler
	 case 1:
		 integrator->euler(dt,in_p2,in_v2);
		break;
	 // semi-implicit euler / symplectic euler
	 case 2:
		 integrator->sympletic_euler(dt, in_p2, in_v2);
		break;
	 // midpoint
	 case 3:
		 integrator->mid_point(dt, in_p2, in_v2);
	 case 4:
	 // semi-implicit backward euler (to be implemented)

		 // TODO
		 break;

	 // analytic
	 case 5:
		 integrator->analytic(dt, in_p2, in_v2);
		 break;
	}


	// set result
	p2 = in_p2;
	v2 = in_v2;
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

// calculate the (specific)force due to a spring of length L between x and the origin
double calculateForce1D(double k, double m, double d, double L, double x, double v){
	return (-k*(L+x) - d*v)/m - g;
}

// calculate the (specific)force due to a spring of length L between x0 and xi
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