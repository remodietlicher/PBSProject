#ifndef SWS_H
#define SWS_H

#define VELOCITY_X 0
#define VELOCITY_Y 1
#define ETA 2

#define INDEX(i, j) (i)*res[0]+(j)

#include <vector>
#include "rigidbody.h"

class SWSolver{
protected:
	float a_ext[2];								// external acceleration
	float v_ext[2];								// external velocity
	float g;									// gravitational constant
	std::vector<float> ground;					// underlying ground
	std::vector<float> eta, height, vel_x, vel_y;
	int res[2];
	float xSize, ySize, dt, dx[2];

public:
	SWSolver(int xRes, int yRes, float xSize, float ySize, float dt);
	const std::vector<float> getHeightMap();
	float getXSize();
	float getYSize();
	int getXRes();
	int getYRes();
	void setEta(std::vector<float> eta);
	void setGround(std::vector<float> ground);
	void setVelocities(std::vector<float> vel_x, std::vector<float> vel_y);
	void setExternalVelocities(float v_ext_x, float v_ext_y);
	void setExternalAccelerations(float a_ext_x, float a_ext_y);
	void advanceTimestep();
	void copyHeights(float*);

protected:
	float interpolate(std::vector<float> &array, float x, float y);
	void calculateHeight();
	void setBoundary();
	void advect(int QUANTITY);
	void updateHeight();
	void updateVelocity();
};

class SWRBSolver : SWSolver {
private:
	Box box;
	std::vector<float> displ_old, displ_new;
public:
	SWRBSolver(int xRes, int yRes, float xSize, float ySize, float dt, Box b);
	void advanceTimestep();
private:
	void handleBodyInteraction();
	std::vector<float> getProjectedIndices();
};


#endif