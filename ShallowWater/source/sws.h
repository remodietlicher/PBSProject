#ifndef SWS_H
#define SWS_H

#define INDEX(i, j) (i)+(j)*grid->xRes

#include <vector>

#include "sw_grid.h"

class SWSolver{
private:
	float a_ext[2];								// external acceleration
	float v_ext[2];								// external velocity
	float g;									// gravitational constant
	float dt, time;
	Sw_grid *grid;		// Geometry

public:
	SWSolver();
	// const std::vector<float> getHeightMap();
	// float getXSize();
	// float getYSize();
	// int getXRes();
	// int getYRes();
	// void setEta(std::vector<float> eta);
	// void setGround(std::vector<float> ground);
	// void setVelocities(std::vector<float> vel_x, std::vector<float> vel_y);
	void setExternalVelocities(float v_ext_x, float v_ext_y);
	void setExternalAccelerations(float a_ext_x, float a_ext_y);
	void advanceTimestep(float _dt);
	// void copyHeights(float*);
	float computeKineticEnergy();
	float computePotentialEnergy();


	void setGrid(Sw_grid *_grid);

private:
	float interpolate(const float *array, float x, float y);

	void calculateHeight();
	void setBoundary();
	void advect(FIELDNAME QUANTITY);
	void updateHeight();
	void updateVelocity();
};

#endif
