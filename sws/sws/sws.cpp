#include "sws.h"

#define INDEX(i, j) i + j*xRes

SWSolver::SWSolver(int xRes, int yRes, float xSize, float ySize, float dt) :
	v_ext({0, 0})
{
	eta.resize(xRes*yRes);
	vel_x.resize(xRes*yRes);
	vel_y.resize(xRes*yRes);
	this->xRes = xRes;
	this->yRes = yRes;
	this->xSize = xSize;
	this->ySize = ySize;
	this->dt = dt;
	this->dx[0] = xSize/xRes;
	this->dx[1] = ySize/yRes;
}

float SWSolver::interpolate(std::vector<float> &quantity, float weight[2], int i0[2]) const {
	double OneMinusWeight[] = {1.0-weight[0], 1.0-weight[1]};
	return (OneMinusWeight[1])*( (OneMinusWeight[0])*quantity[INDEX(i0[0], i0[1])] + weight[0]*quantity[INDEX(i0[0]+1, i0[2])] ) + weight[1]*( (OneMinusWeight[0])*quantity[INDEX(i0[0], i0[1]+1)] + weight[0]*quantity[INDEX(i0[0]+1, i0[1]+1)]);
}

void SWSolver::advect(std::vector<float> array, int QUANTITY){
	float v[2], s, xg[2], xp[2];
	int index;

	std::vector<float> temp;
	temp.resize(xRes*yRes);

	for(int i=1; i<xRes-1; i++)
		for(int j=1; j<yRes-1; j++){
			switch(QUANTITY){
			case 0: // eta (note that this is at the center of a cell)
				v[0] = (vel_x[INDEX(i, j)] + vel_x[INDEX(i+1, j)])*0.5;
				v[1] = (vel_y[INDEX(i, j)] + vel_y[INDEX(i, j+1)])*0.5;
				break;
			case 1: // vel_x (on the left)
				v[0] = vel_x[INDEX(i, j)];
				v[1] = (vel_y[INDEX(i, j)] + vel_y[INDEX(i-1, j)] + vel_y[INDEX(i-1, j+1)] + vel_y[INDEX(i, j+1)])*0.25;
				break;
			case 2: // vel_y (on the bottom)
				v[0] = (vel_x[INDEX(i, j)] + vel_x[INDEX(i-1, j)] + vel_x[INDEX(i-1, j+1)] + vel_x[INDEX(i, j+1)])*0.25;
				v[1] = vel_x[INDEX(i, j)];
				break;
			default: // error
				exit(1);
			}
			// get absolute position of i,j
			float xg[] = {(float)i*dx[0], (float)j*dx[1]};

			// compute biliniear weights
			float xp[2], weight[2]; 					// weight[2]: interpolation weight [0 to 1, where 0 takes x0 and 1 takes xp]
			int i0[2]; 									// index of bottom left gridpoint of the cell, where we interpolate
			for (int i = 0; i < 2; ++i) 				// interpolate in x and y direction
			{
				xp[i] = xg[i] - dt*v[i]; 				// where does the information come from?
				i0[i] = (int)floor(xp[i]/dx[i]);		// index of lower left gridpoint, of the cell of xp
				float x0 = dx[i]*(float)(i0[i]);		// position of i0
				weight[i] = ( xp[i] - x0 ) / dx[i];		// compute weight
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
			eta[INDEX(i, j)] -= eta[INDEX(i, j)] * dt * ((vel_x[INDEX(i+1, j)]-vel_x[INDEX(i, j)])/dx[0] + (vel_y[INDEX(i, j+1)]-vel_y[INDEX(i, j)])/dx[1]);
		}
}

void SWSolver::updateVelocity(){
	
}
