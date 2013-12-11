#include "sws.h"
#include "vec3f.h"
#include <iostream>
#include <math.h>


SWSolver::SWSolver() {
	// initialize defaults
	a_ext[0] = 0;
	a_ext[1] = 0;
	v_ext[0] = 0;
	v_ext[1] = 0;
	g = 9.81f;
	time = 0.0;

	// calculateHeight();
}

void SWSolver::setGrid(Sw_grid *_grid){
	grid = _grid;

	calculateHeight();
}

float SWSolver::interpolate(const float *array, float x, float y) {
	const int I = (int)x;
	const int J = (int)y;
	const float s1 = x - I;
	const float s0 = 1.f - s1;
	const float t1 = y - J;
	const float t0 = 1.f-t1;
	return  s0*(t0*array[INDEX(I, J)]+t1*array[INDEX(I, J+1)] )+s1*(t0*array[INDEX(I+1, J)] +t1*array[INDEX(I+1, J+1)]);
}

void SWSolver::advect(FIELDNAME QUANTITY){
	float v[2];

	const float *array = grid->oldFields[QUANTITY];
	const float *vel_x = grid->oldFields[VELX];
	const float *vel_y = grid->oldFields[VELY];

	float *noi  = grid->newFields[QUANTITY];

	float dx = grid->dx;
	int xRes = grid->xRes;
	int yRes = grid->yRes;

	for(int i=1; i<xRes-1; i++)
		for(int j=1; j<yRes-1; j++){
			// set external velocities present at all times
			v[0] = v_ext[0];
			v[1] = v_ext[1];
			switch(QUANTITY){
			case ETA: // eta (note that this is at the center of a cell)
				v[0] += (vel_x[INDEX(i, j)] + vel_x[INDEX(i+1, j)])*0.5f;
				v[1] += (vel_y[INDEX(i, j)] + vel_y[INDEX(i, j+1)])*0.5f;
				break;
			case VELX: // vel_x (on the left)
				v[0] += vel_x[INDEX(i, j)];
				v[1] += (vel_y[INDEX(i, j)] + vel_y[INDEX(i-1, j)] + vel_y[INDEX(i-1, j+1)] + vel_y[INDEX(i, j+1)])*0.25f;
				break;
			case VELY: // vel_y (on the bottom)
				v[0] += (vel_x[INDEX(i, j)] + vel_x[INDEX(i-1, j)] + vel_x[INDEX(i-1, j+1)] + vel_x[INDEX(i, j+1)])*0.25f;
				v[1] += vel_y[INDEX(i, j)];
				break;
			default: // error
				printf("Error in SWSolver::advect\n\twrong FIELD=%d",QUANTITY);
				exit(1);
			}

			// backtrace position
			float srcpi = (float)i - v[0] * dt * 1.0/dx;
			float srcpj = (float)j - v[1] * dt * 1.0/dx;

			// clamp range of accesses (dont access boundaries (at 0 and res-1))
			if(srcpi<0.) srcpi = 1.f;
			if(srcpj<0.) srcpj = 1.f;
			if(srcpi>xRes-1.) srcpi = xRes-2.f;
			if(srcpj>yRes-1.) srcpj = yRes-2.f;

			// interpolate source value
			noi[INDEX(i, j)] = interpolate(array, srcpi, srcpj);
		}
		// grid->switchOldNew(QUANTITY);
}

void SWSolver::updateHeight(){
	int xRes = grid->xRes;
	int yRes = grid->yRes;
	float dx = grid->dx;

	float *eta    = grid->oldFields[ETA];
	float *height = grid->oldFields[HEIGHT];
	float *ground  = grid->oldFields[GROUND];
	const float *vel_x = grid->oldFields[VELX];
	const float *vel_y = grid->oldFields[VELY];

	for (int i = 1; i < xRes - 1; i++)
	for (int j = 1; j < yRes - 1; j++){	
		eta[INDEX(i, j)] += -eta[INDEX(i, j)] * dt * ((vel_x[INDEX(i + 1, j)] - vel_x[INDEX(i, j)]) / dx + (vel_y[INDEX(i, j + 1)] - vel_y[INDEX(i, j)]) / dx);
		if (eta[INDEX(i, j)] < 0.0) 
			eta[INDEX(i, j)] = 0.0;
		height[INDEX(i, j)] = eta[INDEX(i, j)] + ground[INDEX(i, j)];
	}
	
}

void SWSolver::updateVelocity(){
	int xRes = grid->xRes;
	int yRes = grid->yRes;
	float dx = grid->dx;
	float *vel_x = grid->oldFields[VELX];
	float *vel_y = grid->oldFields[VELY];
	float *height = grid->oldFields[HEIGHT];

	for(int i=2; i<xRes-1; i++)
		for(int j=1; j<yRes-1; j++){
			vel_x[INDEX(i, j)] += (g*((height[INDEX(i-1, j)]-height[INDEX(i, j)])/dx)+a_ext[0])*dt;
		}
	for(int i=1; i<xRes-1; i++)
		for(int j=2; j<yRes-1; j++){
			vel_y[INDEX(i, j)] += (g*((height[INDEX(i, j-1)]-height[INDEX(i, j)])/dx)+a_ext[0])*dt;
		}
}

void SWSolver::setBoundary(){
	int xRes = grid->xRes;
	int yRes = grid->yRes;
	float *height = grid->oldFields[HEIGHT];

	for(int i=0; i<xRes; i++){
		height[INDEX(i, 0)]      = height[INDEX(i, 1)];				// lower boundary
		height[INDEX(i, yRes-1)] = height[INDEX(i, yRes-2)];	// upper boundary
	}
	for(int j=0; j<yRes; j++){
		height[INDEX(0, j)]      = height[INDEX(1, j)];				// left boundary
		height[INDEX(xRes-1, j)] = height[INDEX(xRes-2, j)];	// right boundary
	}
}

void SWSolver::advanceTimestep(float time_goal){
		float *eta = grid->oldFields[ETA];
		int xRes = grid->xRes;
		int yRes = grid->yRes;
		float H = 0.;
		for (int i = 0; i < xRes - 1; i++)
		for (int j = 0; j < yRes - 1; j++)
			H += eta[i+j*xRes];
		H /= (float)(xRes*yRes);
		float maxdt = 0.1*grid->dx / sqrt(g*H);
		
		while(time < time_goal){
			
			dt = std::min(maxdt, time_goal - time);

			advect(ETA);
			advect(VELX);
			advect(VELY);
			grid->switchOldNew(ETA);
			grid->switchOldNew(VELX);
			grid->switchOldNew(VELY);
			updateHeight();
			updateVelocity();
			setBoundary();
			time += dt;
		}

}

void SWSolver::setExternalVelocities(float v_ext_x, float v_ext_y){
	this->v_ext[0] = v_ext_x;
	this->v_ext[1] = v_ext_y;
}
void SWSolver::setExternalAccelerations(float a_ext_x, float a_ext_y){
	this->a_ext[0] = a_ext_x;
	this->a_ext[1] = a_ext_y;
}
void SWSolver::calculateHeight(){
	float *height = grid->oldFields[HEIGHT];
	float *eta = grid->oldFields[ETA];
	float *ground = grid->oldFields[GROUND];
	float xRes = grid->xRes;
	float yRes = grid->yRes;
	for(int i=0; i<xRes; i++)
		for(int j=0; j<yRes; j++){
			height[INDEX(i, j)] = eta[INDEX(i, j)] + ground[INDEX(i, j)];
		}
}

float SWSolver::computeKineticEnergy(){
	float *vel_x = grid->oldFields[VELX];
	float *vel_y = grid->oldFields[VELY];
	float energy = 0.;
	for (int i = 0; i < grid->xRes; ++i)
	for (int j = 0; j < grid->yRes; ++j)
	{
		float u = vel_x[i+j*grid->xRes];
		float v = vel_y[i+j*grid->xRes];
		energy += 0.5 * (u*u + v*v);
	}
	return energy;
}

float SWSolver::computePotentialEnergy(){
	float *height = grid->oldFields[HEIGHT];
	float energy = 0.;
	for (int i = 0; i < grid->xRes; ++i)
	for (int j = 0; j < grid->yRes; ++j){
		energy += g*height[i+j*grid->xRes];
	}
	return energy;
}
