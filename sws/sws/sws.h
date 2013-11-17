#ifndef SWS_H
#define SWS_H

#define VELOCITY_X 0
#define VELOCITY_Y 1
#define ETA 2

#include <vector>

class SWSolver{
private:
	const float a_ext;
	const float v_ext[2];
	std::vector<float> eta, vel_x, vel_y;
	int xRes, yRes;
	float xSize, ySize, dt, dx[2];
public:
	SWSolver(int xRes, int yRes, float xSize, float ySize, float dt);
	void advect(std::vector<float> array, int QUANTITY);
	void updateHeight();
	void updateVelocity();
private:
	float interpolate(std::vector<float> &quantity, float weight[2], int i0[2]) const;
};

#endif