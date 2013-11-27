#include "sws.h"
#include "vec3f.h"
#include <iostream>


SWSolver::SWSolver(int xRes, int yRes, float xSize, float ySize, float dt) {
	// initialize defaults
	a_ext[0] = 0;
	a_ext[1] = 0;
	v_ext[0] = 0;
	v_ext[1] = 0;
	g = 9.81f;
	ground.resize(xRes*yRes);	// initializes all values to zero
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

	std::vector<float> temp, *array;
	temp.resize(res[0]*res[1]);

	for(int i=1; i<res[0]-1; i++)
		for(int j=1; j<res[1]-1; j++){
			// set external velocities present at all times
			v[0] = v_ext[0];
			v[1] = v_ext[1];
			switch(QUANTITY){
			case 0: // eta (note that this is at the center of a cell)
				array = &eta;
				v[0] += (vel_x[INDEX(i, j)] + vel_x[INDEX(i+1, j)])*0.5f;
				v[1] += (vel_y[INDEX(i, j)] + vel_y[INDEX(i, j+1)])*0.5f;
				break;
			case 1: // vel_x (on the left)
				array = &vel_x;
				v[0] += vel_x[INDEX(i, j)];
				v[1] += (vel_y[INDEX(i, j)] + vel_y[INDEX(i-1, j)] + vel_y[INDEX(i-1, j+1)] + vel_y[INDEX(i, j+1)])*0.25f;
				break;
			case 2: // vel_y (on the bottom)
				array = &vel_y;
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
			temp[INDEX(i, j)] = interpolate(*array, srcpi, srcpj);
		}
		for (int i = 0; i < res[0]; i++)
			for (int j = 0; j < res[1]; j++){
				(*array)[INDEX(i, j)] = temp[INDEX(i, j)];
	}
}

void SWSolver::updateHeight(){
	//std::cout << "update height:" << std::endl;
	for (int i = 1; i < res[0] - 1; i++) {
		//std::cout << std::endl << i << ": ";
		for (int j = 1; j < res[1] - 1; j++){

			float oldEta = eta[INDEX(i, j)];

			eta[INDEX(i, j)] += -eta[INDEX(i, j)] * dt * ((vel_x[INDEX(i + 1, j)] - vel_x[INDEX(i, j)]) / dx[0] + (vel_y[INDEX(i, j + 1)] - vel_y[INDEX(i, j)]) / dx[1]);
			if (eta[INDEX(i, j)] < ground[INDEX(i, j)]) eta[INDEX(i, j)] = ground[INDEX(i, j)];
			height[INDEX(i, j)] = eta[INDEX(i, j)] + ground[INDEX(i, j)];

			float newEta = eta[INDEX(i, j)];

			//std::cout << newEta - oldEta<<" ";
		}
		
	}
	//std::cout <<std::endl;
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
		height[INDEX(i, 0)] = height[INDEX(i, 1)];				// lower boundary
//		height[INDEX(i, 0)] = 7.0;				// lower boundary
		height[INDEX(i, res[0]-1)] = height[INDEX(i, res[0]-2)];	// upper boundary
	}
	for(int j=0; j<res[1]; j++){
		height[INDEX(0, j)] = height[INDEX(1, j)];				// left boundary
//		height[INDEX(0, j)] = 3.0;				// left boundary
		height[INDEX(res[1]-1, j)] = height[INDEX(res[1]-2, j)];	// right boundary
	}
}

void SWSolver::advanceTimestep(){
		//std::cout << "advecting eta..." << std::endl;
		advect(ETA);
		//std::cout << "advecting velocity_x..." << std::endl;
		advect(VELOCITY_X);
		//std::cout << "advecting velocity_y..." << std::endl;
		advect(VELOCITY_Y);

		//std::cout << "updating heights..." << std::endl;
		updateHeight();

		//std::cout << "updating velocities..." << std::endl;
		updateVelocity();

		//std::cout << "setting boundaries..." << std::endl;
		setBoundary();
}

void SWSolver::copyHeights(float* buffer)
{
	for (int i = 1; i < res[0] - 1; i++) {
		//std::cout << std::endl << i << ": ";
		for (int j = 1; j < res[1] - 1; j++){

			float h = height[INDEX(i, j)];
			buffer[INDEX(i, j)] = h;
		}
	}
}

const std::vector<float> SWSolver::getHeightMap(){
	return height;
}

void SWSolver::setEta(std::vector<float> eta){
	this->eta = eta;
	calculateHeight();
}
void SWSolver::setGround(std::vector<float> ground){
	this->ground = ground;
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
			height[INDEX(i, j)] = eta[INDEX(i, j)] + ground[INDEX(i, j)];
		}
}

float SWSolver::getXSize(){
	return xSize;
}
float SWSolver::getYSize(){
	return ySize;
}
int SWSolver::getXRes(){
	return res[0];
}
int SWSolver::getYRes(){
	return res[1];
}

SWRBSolver::SWRBSolver(int xRes, int yRes, float xSize, float ySize, float dt, Box b) : 
	SWSolver(xRes, yRes, xSize, ySize, dt),
	box(b)
{
	displ_old.resize(xRes*yRes);
	displ_new.resize(xRes*yRes);
}

void SWRBSolver::advanceTimestep(){
	std::vector<float> projIndices = getProjectedIndices();
	for(int i=0; i<projIndices.size(); i++){
		
	}
	SWSolver::advanceTimestep();

	cout << "handling body interaction..." << endl;

	handleBodyInteraction();
}

std::vector<float> SWRBSolver::getProjectedIndices(){
	std::vector<float> indices;
	indices.resize(1);
	return indices;
}

void SWRBSolver::handleBodyInteraction(){
	
}
