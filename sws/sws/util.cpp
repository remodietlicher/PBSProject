#include "util.h"
#include "assert.h"

Grid2Df::Grid2Df(int Nx, int Ny){
	this->Nx = Nx;
	this->Ny = Ny;
	grid.resize(Nx*Ny);
}

float& Grid2Df::operator()(int i, int j){
	assert(i<Nx && i>=0 && j<Ny && j>=0);
	return grid[i*Nx + j];
}
