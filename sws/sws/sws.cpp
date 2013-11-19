#include "sws.h"
#include <iostream>

SWSolver::SWSolver(int xRes, int yRes, float xSize, float ySize, float dt) {
	// initialize defaults
	a_ext[0] = 0;
	a_ext[1] = 0;
	v_ext[0] = 0;
	v_ext[1] = 0;
	g = 0.81f;
	terrain.resize(xRes*yRes);	// initializes all values to zero
	eta.resize(xRes*yRes);
	eta.assign(xRes*yRes, 1);	// put some default values
	height.resize(xRes*yRes);
	vel_x.resize(xRes*yRes);
	vel_x.assign(xRes*yRes, 0.2f);
	vel_y.resize(xRes*yRes);
	// rest
	this->xRes = xRes;
	this->yRes = yRes;
	this->xSize = xSize;
	this->ySize = ySize;
	this->dt = dt;
	this->dx[0] = xSize/xRes;
	this->dx[1] = ySize/yRes;

	calculateHeight();
}

float SWSolver::interpolate(std::vector<float> &quantity, float weight[2], int i0[2]) const {
	float OneMinusWeight[] = {1.0f-weight[0], 1.0f-weight[1]};
	return (OneMinusWeight[1])*( (OneMinusWeight[0])*quantity[INDEX(i0[0], i0[1])] + weight[0]*quantity[INDEX(i0[0]+1, i0[1])] ) + weight[1]*( (OneMinusWeight[0])*quantity[INDEX(i0[0], i0[1]+1)] + weight[0]*quantity[INDEX(i0[0]+1, i0[1]+1)]);
}

void SWSolver::advect(int QUANTITY){
	float v[2];

	std::vector<float> temp, array;
	temp.resize(xRes*yRes);

	for(int i=1; i<xRes-1; i++)
		for(int j=1; j<yRes-1; j++){
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
			// get absolute position of i,j
			float xg[] = {(float)i*dx[0], (float)j*dx[1]};

			// compute biliniear weights
			float xp[2], weight[2]; 					// weight[2]: interpolation weight [0 to 1, where 0 takes x0 and 1 takes xp]
			int i0[2]; 									// index of bottom left gridpoint of the cell, where we interpolate
			for (int l = 0; l < 2; l++) 				// interpolate in x and y direction
			{
				xp[l] = xg[l] - dt*v[l]; 				// where does the information come from?
				i0[l] = (int)floor(xp[l]/dx[l]);		// index of lower left gridpoint, of the cell of xp
				float x0 = dx[l]*(float)(i0[l]);		// position of i0
				weight[l] = ( xp[l] - x0 ) / dx[l];		// compute weight
			}

			temp[INDEX(i, j)] = interpolate(array, weight, i0);	
		}
		for (int y = 0; y < yRes; y++)
			for (int x = 0; x < xRes; x++){
				array[x + y*xRes] = temp[x + y*xRes];
	}
}

void SWSolver::updateHeight(){
	for(int i=1; i<xRes-1; i++)
		for(int j=1; j<yRes-1; j++){
			height[INDEX(i, j)] = eta[INDEX(i, j)] + terrain[INDEX(i, j)];
			height[INDEX(i, j)] -= height[INDEX(i, j)] * dt * ((vel_x[INDEX(i+1, j)]-vel_x[INDEX(i, j)])/dx[0] + (vel_y[INDEX(i, j+1)]-vel_y[INDEX(i, j)])/dx[1]);
		}
}

void SWSolver::updateVelocity(){
	for(int i=1; i<xRes-1; i++)
		for(int j=1; j<yRes-1; j++){
			if(i>1) vel_x[INDEX(i, j)] = (-g*((height[INDEX(i-1, j)]-height[INDEX(i, j)])/dx[0])+a_ext[0])*dt;
			if(j>1) vel_y[INDEX(i, j)] = (-g*((height[INDEX(i, j-1)]-height[INDEX(i, j)])/dx[0])+a_ext[0])*dt;
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
	for(int i=0; i<xRes; i++)
		for(int j=0; j<yRes; j++){
			height[INDEX(i, j)] = eta[INDEX(i, j)] + terrain[INDEX(i, j)];
		}
}
