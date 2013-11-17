#ifndef SWS_H
#define SWS_H

#define VELOCITY_X 0
#define VELOCITY_Y 1
#define ETA 2

#include <vector>

class SWSolver{
private:
	std::vector<float> eta, vel_x, vel_y;
	int xRes, yRes;
	float xSize, ySize, dt, dx[2];
public:
	SWSolver(int xRes, int yRes, float xSize, float ySize, float dt);
	void advect(int QUANTITY);
};

#endif