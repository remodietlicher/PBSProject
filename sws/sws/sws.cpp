#include "sws.h"
#include <iostream>

SWSolver::SWSolver(int xRes, int yRes, float xSize, float ySize, float dt) {
	// initialize defaults
	a_ext[0] = 0;
	a_ext[1] = 0;
	v_ext[0] = 0;
	v_ext[1] = 0;
	g = -9.81f;
	terrain.resize(xRes*yRes);	// initializes all values to zero
	eta.resize(xRes*yRes);
	eta.assign(xRes*yRes, 1);	// put some default values
	height.resize(xRes*yRes);
	vel_x.resize(xRes*yRes);
//	vel_x.assign(xRes*yRes, 5.0f);
	vel_y.resize(xRes*yRes);
//	vel_x.assign(xRes*yRes, 5.0f);
	// rest
	this->res[0] = xRes;
	this->res[1] = yRes;
	this->xSize = xSize;
	this->ySize = ySize;
	this->dt = dt;
	this->dx[0] = xSize/xRes;
	this->dx[1] = ySize/yRes;

	calculateHeight();
}

float SWSolver::interpolate(std::vector<float> &array, float x, float y) {
	const int I = (int)x;
	const int J = (int)y;
	const float s1 = x - I;
	const float s0 = 1.f - s1;
	const float t1 = y - J;
	const float t0 = 1.f-t1;
	return  s0*(t0*array[INDEX(I, J)]+t1*array[INDEX(I, J+1)] )+s1*(t0*array[INDEX(I+1, J)] +t1*array[INDEX(I+1, J+1)]);
}

void SWSolver::advect(int QUANTITY){
	float v[2];

	std::vector<float> temp, array;
	temp.resize(res[0]*res[1]);

	for(int i=1; i<res[0]-1; i++)
		for(int j=1; j<res[1]-1; j++){
			// set external velocities present at all times
			v[0] = v_ext[0];
			v[1] = v_ext[1];
			switch(QUANTITY){
			case 0: // eta (note that this is at the center of a cell)
				array = eta;
				v[0] += (vel_x[INDEX(i, j)] + vel_x[INDEX(i+1, j)])*0.5f;
				v[1] += (vel_y[INDEX(i, j)] + vel_y[INDEX(i, j+1)])*0.5f;
				break;
			case 1: // vel_x (on the left)
				array = vel_x;
				v[0] += vel_x[INDEX(i, j)];
				v[1] += (vel_y[INDEX(i, j)] + vel_y[INDEX(i-1, j)] + vel_y[INDEX(i-1, j+1)] + vel_y[INDEX(i, j+1)])*0.25f;
				break;
			case 2: // vel_y (on the bottom)
				array = vel_y;
				v[0] += (vel_x[INDEX(i, j)] + vel_x[INDEX(i-1, j)] + vel_x[INDEX(i-1, j+1)] + vel_x[INDEX(i, j+1)])*0.25f;
				v[1] += vel_x[INDEX(i, j)];
				break;
			default: // error
				exit(1);
			}

			// backtrace position
			float srcpi = (float)i - v[0] * dt * 1/dx[0];
			float srcpj = (float)j - v[1] * dt * 1/dx[1];

			// clamp range of accesses (dont access boundaries (at 0 and res-1))
			if(srcpi<0.) srcpi = 1.f;
			if(srcpj<0.) srcpj = 1.f;
			if(srcpi>res[0]-1.) srcpi = res[0]-2.f;
			if(srcpj>res[1]-1.) srcpj = res[1]-2.f;

			// interpolate source value
			temp[INDEX(i, j)] = interpolate(array, srcpi, srcpj);
/*			// get absolute position of i,j
			float xg[] = {(float)i*dx[0], (float)j*dx[1]};

			// compute biliniear weights
			float xp[2], weight[2], x0[2];				// weight[2]: interpolation weight [0 to 1, where 0 takes x0 and 1 takes xp]
			int i0[2]; 									// index of bottom left gridpoint of the cell, where we interpolate
			for (int l = 0; l < 2; l++) 				// interpolate in x and y direction
			{
				xp[l] = xg[l] - dt*v[l]; 				// where does the information come from?
				i0[l] = (int)floor(xp[l]/dx[l]);		// index of lower left gridpoint, of the cell of xp
				if(i0[l]<0.0f) i0[l] = 0.0f;			// clamp at border
				if(i0[l]>res[l]-1) i0[l] = (float)(res[l]-1);
				x0[l] = dx[l]*(float)(i0[l]);			// position of i0
				weight[l] = ( xp[l] - x0[l] ) / dx[l];	// compute weight
			}

			temp[INDEX(i, j)] = interpolate(array, weight, i0);	*/
		}
		for (int i = 0; i < res[0]; i++)
			for (int j = 0; j < res[1]; j++){
				array[INDEX(i, j)] = temp[INDEX(i, j)];
	}
}

void SWSolver::updateHeight(){
	for(int i=1; i<res[0]-1; i++)
		for(int j=1; j<res[1]-1; j++){
			eta[INDEX(i, j)] += -eta[INDEX(i, j)] * dt * ((vel_x[INDEX(i+1, j)]-vel_x[INDEX(i, j)])/dx[0] + (vel_y[INDEX(i, j+1)]-vel_y[INDEX(i, j)])/dx[1]);
			height[INDEX(i, j)] = eta[INDEX(i, j)] + terrain[INDEX(i, j)];
		}
}

void SWSolver::updateVelocity(){
	for(int i=2; i<res[0]-1; i++)
		for(int j=1; j<res[1]-1; j++){
			vel_x[INDEX(i, j)] += (g*((height[INDEX(i-1, j)]-height[INDEX(i, j)])/dx[0])+a_ext[0])*dt;
		}
	for(int i=1; i<res[0]-1; i++)
		for(int j=2; j<res[1]-1; j++){
			vel_y[INDEX(i, j)] += (g*((height[INDEX(i, j-1)]-height[INDEX(i, j)])/dx[0])+a_ext[0])*dt;
		}
}

void SWSolver::setBoundary(){
	for(int i=0; i<res[0]; i++){
		eta[INDEX(i, 0)] = eta[INDEX(i, 1)];				// lower boundary
		eta[INDEX(i, res[0]-1)] = eta[INDEX(i, res[0]-2)];	// upper boundary
	}
	for(int j=0; j<res[1]; j++){
		eta[INDEX(0, j)] = eta[INDEX(1, j)];				// left boundary
		eta[INDEX(res[1]-1, j)] = eta[INDEX(res[1]-2, j)];	// right boundary
	}
}

std::vector<float> SWSolver::getHeightMap(){
	return height;
}

void SWSolver::setEta(std::vector<float> eta){
	this->eta = eta;
	calculateHeight();
}
void SWSolver::setTerrain(std::vector<float> terrain){
	this->terrain = terrain;
	calculateHeight();
}
void SWSolver::setVelocities(std::vector<float> vel_x, std::vector<float> vel_y){
	this->vel_x = vel_x;
	this->vel_y = vel_y;
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
	for(int i=0; i<res[0]; i++)
		for(int j=0; j<res[1]; j++){
			height[INDEX(i, j)] = eta[INDEX(i, j)] + terrain[INDEX(i, j)];
		}
}
