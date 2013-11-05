#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include "solver.h"

class SolverGaussSeidel : public SolverBase 
{
	public:
		SolverGaussSeidel(int iterations, double dt, double acc) :
//SolverGaussSeidel::SolverGaussSeidel(int iterations, double dt, double acc) :
  SolverBase(iterations, dt, acc) { };

		~SolverGaussSeidel() {};

		void solvePoisson(int xRes, int yRes, double* field, double* b);

	private:
};

#endif
