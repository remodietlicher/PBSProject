#ifndef SWS_H
#define SWS_H

#define VELOCITY_X 0
#define VELOCITY_Y 1
#define ETA 2

#define INDEX(i, j) (i)+(j)*xRes
#include <vector>

class SWSolver{
private:
	float a_ext[2];								// external acceleration
	float v_ext[2];								// external velocity
	float g;									// gravitational constant
	std::vector<float> terrain;					// underlying terrain
	std::vector<float> eta, height, vel_x, vel_y;
	int xRes, yRes;
	float xSize, ySize, dt, dx[2];
public:
	SWSolver(int xRes, int yRes, float xSize, float ySize, float dt);
	void advect(int QUANTITY);
	void updateHeight();
	void updateVelocity();
	std::vector<float> getHeightMap();
	void setEta(std::vector<float> eta);
	void setTerrain(std::vector<float> terrain);
	void setVelocities(std::vector<float> vel_x, std::vector<float> vel_y);
	void setExternalVelocities(float v_ext_x, float v_ext_y);
	void setExternalAccelerations(float a_ext_x, float a_ext_y);
private:
	float interpolate(std::vector<float> &quantity, float weight[2], int i0[2]) const;
	void calculateHeight();
};

#endif