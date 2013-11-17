#ifndef UTIL_H
#define UTIL_H

#include <vector>

class Grid2Df{
private:
	int Nx, Ny;
	std::vector<float> grid;
public:
	Grid2Df();
	Grid2Df(int Nx, int Ny);
	float &operator()(int i, int j);
};

#endif