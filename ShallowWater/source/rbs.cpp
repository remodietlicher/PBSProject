#include "rbs.h"

RigidBodySolver::RigidBodySolver(Box *b) :
	body(b)
	{}

Matrix3x3f RigidBodySolver::star(Vector3f v){
	return 	Matrix3x3f(0.0f, -v[2], v[1],
                       v[2], 0.0f, -v[0],
                       -v[1], v[0], 0.0f);
}

//! parameter: dt: timestep, F: array of external forces, r: where the external forces act on body, N: number of ext. F.
void RigidBodySolver::advanceTimestep(float dt, std::vector<Vector3f> F, std::vector<Vector3f> r){
	sumExternalForces(F, r);
	performEulerIntegration(dt);
	setBodyBox();
}

//! parameter: r[i] is the position on the body where F[i] acts
void RigidBodySolver::sumExternalForces(std::vector<Vector3f> F, std::vector<Vector3f> r){
	body->force = Vector3f(0.0f, 0.0f, 0.0f);
	body->torque = Vector3f(0.0f, 0.0f, 0.0f);
	int N = (int)F.size();
	for(int i=0; i<N; i++){
		body->force += F[i];
        
        Vector3f torque = (r[i] - body->x).cross(F[i]);
            body->torque += torque;
//        if(body->torque.length() > 0.01)
//            body->torque = body->torque/body->torque.length() * 0.01;
            //        printf("torque= %f, F= %f, r= %f\n", torque.length(), F[i].length(), (r[i] - body->x).length());
	}
}

void RigidBodySolver::performEulerIntegration(float dt){
	body->x    += dt*body->v;
	body->v    += 1/body->mass*dt*body->force;
	body->R    += dt*star(body->omega)*body->R;
	body->Iinv  = body->R * body->Ibodyinv * body->R.transposed();
	body->L    += dt*body->torque;
	body->omega = body->Iinv*body->L;
    if(body->L.length() > 0.1)
        body->L = body->L/body->L.length() * 0.1;

}

void RigidBodySolver::performImplicitEulerIntegration(float dt){
	// TODO adjust for implicit Euler integration
	body->v += 1/body->mass*dt*body->force;
	body->x += dt*body->v;
	body->R += dt*star(body->omega)*body->R;
	body->Iinv = body->R*body->Ibodyinv*body->R.transposed();
	body->L += dt*body->torque;
	body->omega = body->Iinv*body->L;
}

void RigidBodySolver::setBodyBox(){
	body->x0 = body->R*body->x0;
	body->y0 = body->R*body->y0;
	body->z0 = body->R*body->z0;
}

Box* RigidBodySolver::getBody(){
	return body;
}