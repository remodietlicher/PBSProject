#ifndef SWS_H
#define SWS_H

#define INDEX(i, j) (i)+(j)*grid->xRes

#include <vector>

#include "sw_grid.h"
#include "rigidbody.h"
#include "rbs.h"

class SWSolver{
protected:
	float a_ext[2];								// external acceleration
	float v_ext[2];								// external velocity
	float dt, time;
	Sw_grid *grid;		// Geometry

public:
	SWSolver();
	void setExternalVelocities(float v_ext_x, float v_ext_y);
	void setExternalAccelerations(float a_ext_x, float a_ext_y);
	void advanceTimestep(float _dt);
	float computeKineticEnergy();
	float computePotentialEnergy();
	void setGrid(Sw_grid *_grid);
	float g;									// gravitational constant

protected:
	float interpolate(const float *array, float x, float y);

	void calculateHeight();
	void setBoundary();
	void advect(FIELDNAME QUANTITY);
	void updateHeight();
	void updateVelocity();
};

class SWRBSolver {
private:
	Box *box;
	RigidBodySolver rbs;
	std::vector<float> displ_old, displ_new;
	float alpha, rho, g;
	Sw_grid *grid;

public:
	SWRBSolver(Sw_grid *_grid, Box *_box);
	void advanceTimestep(float dt);
	Box* getBody();
	std::vector<float>* getDisplacement();
	void testSorting(); // only for debug purposes
private:
	void handleBodyInteraction();
	void setDisplacement(std::vector<Vector3f> &positions, int x_min, int x_max, int y_min, int y_max, bool* isIntersecting);
	void estimateIndices(Vector3f vertices[8], int &x_min, int &x_max, int &y_min, int &y_max);
	std::vector<Vector3f> getConvexHull8XY(Vector3f *vertices);
	void bubbleSortVert(int coord, Vector3f *A, int n);
	bool calculateDisplacement(int i, int j, float &b, Vector3f &r);
};

#endif
