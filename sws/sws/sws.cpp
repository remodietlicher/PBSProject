#include "sws.h"

#define INDEX(i, j) i + j*yRes

SWSolver::SWSolver(int xRes, int yRes, float xSize, float ySize, float dt){
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

void SWSolver::advect(int QUANTITY){
	float v[2], s, xg[2], xp[2];
	int index;
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

			// update with interpolation weights
			int x0y0 = i0[0]   + ( i0[1]  )*xRes; // lower left
			int x1y0 = i0[0]+1 + ( i0[1]  )*xRes; // lower right
			int x0y1 = i0[0]   + ( i0[1]+1)*xRes; // upper left
			int x1y1 = i0[0]+1 + ( i0[1]+1)*xRes; // upper right
			double OneMinusWeight[] = {1.0-weight[0], 1.0-weight[1]};
			xVelocityTemp[x + y*xRes] = (OneMinusWeight[1])*( (OneMinusWeight[0])*xVelocity[x0y0] + weight[0]*xVelocity[x1y0] ) + weight[1]*( (OneMinusWeight[0])*xVelocity[x0y1] + weight[0]*xVelocity[x1y1] );
			yVelocityTemp[x + y*xRes] = (OneMinusWeight[1])*( (OneMinusWeight[0])*yVelocity[x0y0] + weight[0]*yVelocity[x1y0] ) + weight[1]*( (OneMinusWeight[0])*yVelocity[x0y1] + weight[0]*yVelocity[x1y1] );
			densityTemp[x + y*xRes]   = (OneMinusWeight[1])*( (OneMinusWeight[0])*density[x0y0]   + weight[0]*density[x1y0] )   + weight[1]*( (OneMinusWeight[0])*density[x0y1]   + weight[0]*density[x1y1] );
		}
}
