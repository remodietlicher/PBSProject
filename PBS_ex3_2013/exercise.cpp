#include "gauss_seidel.h"
#include "fluid2d.h"
#include "util.h"

#define GETF(x,y) field[(x) + (y) * xRes]
#define GETB(x,y) b[(x) + (y) * xRes]


// Problem 1
#define INDEX(x,y) ((x) + (y) * xRes)
// size of the Poisson finite differences matrix
#define POM_SIZE (SIM_RES * SIM_RES)
#define INIT_VECTOR(vec) for(int _i=0;_i<POM_SIZE;_i++) {vec[_i] = 0.0;} 

// _pressure[x]  <pressure for cell x>
// _divergence[x] <divergence for cell x>
void ExSolvePoisson(int xRes, int yRes, int _iterations, double _accuracy, double* _pressure, double* _divergence)
{
	static int matrix[POM_SIZE][POM_SIZE];
	static bool init = false;
	static double dxx = pow(1.0 / SIM_RES, 2);

	if (!init) {
		// init to 0
		for (int n = 0; n < POM_SIZE; n++) {
			for (int m = 0; m < POM_SIZE; m++) {
				matrix[m][n] = 0;
			}
		}

		for (int y = 1; y < yRes - 1; y++) {
			for (int x = 1; x < xRes - 1; x++)
			{
				int idx = INDEX(x, y);//index of cell x,y
				matrix[idx][idx] = 4;

				if (x + 1 < xRes - 1) {
					matrix[idx][INDEX(x+1, y)] = -1;
				}
				if (y + 1 < yRes - 1) {
					matrix[idx][INDEX(x, y+1)] = -1;
				}
				if (x - 1 > 0) {
					matrix[idx][INDEX(x-1, y)] = -1;
				}
				if (y - 1 > 0) {
					matrix[idx][INDEX(x, y-1)] = -1;
				}
			}
		}
		init = true;
	}

	double* p_old = new double[POM_SIZE];
	double* p_new = new double[POM_SIZE];
	double* div = new double[POM_SIZE];
	INIT_VECTOR(p_old);
	INIT_VECTOR(p_new);
	INIT_VECTOR(div);

	// note that the boundaries are handles by the framework, so you iterations should be similar to:
	for (int y = 1; y < yRes - 1; y++) {
		for (int x = 1; x < xRes - 1; x++)
		{
			p_old[INDEX(x, y)] = _pressure[INDEX(x, y)];
			//GETF(x,y) ... GETB(x,y)
		}
	}

	for (int iter = 0; iter < _iterations;iter++) {

		for (int k = 0; k < POM_SIZE;k++) {
			double sum = 0.0;
			for (int i = 0; i < POM_SIZE; i++) {

				if (i < k) {
					sum += matrix[k][i] * p_new[i];
				}
				else if (i>k){
					sum += matrix[k][i] * p_old[i];
				}
			}
			p_new[k] = dxx / matrix[k][k] *(div[k] - sum);
		}

		// get residual
		double residual = 0.0;
		for (int k = 0; k < POM_SIZE; k++) {
			residual += pow(p_new[k] - p_old[k],2);
		}
		residual = sqrt(residual);
		if (iter == _iterations - 1) {
			printf("Pressure solver: iter=%d , res=%f \n", iter, residual);
		}
		if (residual < _accuracy) {
			printf("Pressure solver: iter=%d , converged \n", iter, residual);
			break; 
		}
	}

	//update pressur
	for (int y = 1; y < yRes - 1; y++) {
		for (int x = 1; x < xRes - 1; x++)
		{
			_pressure[INDEX(x, y)] = p_new[INDEX(x, y)];
			//GETF(x,y) ... GETB(x,y)
		}
	}

	//cleanup
	delete p_new;
	delete p_old;
	delete div;

	// for your debugging, and ours, please add these prints after every iteration
	/*if(iter==_iterations-1) 
		printf("Pressure solver: iter=%d , res=%f \n",iter, residual);
	if(residual<_accuracy) {
		printf("Pressure solver: iter=%d , converged \n",iter,residual);
		break; // optional
	} */

}

// Probelm 2
void ExCorrectVelocities(int _xRes, int _yRes, double _dt, const double* _pressure, double* _xVelocity, double* _yVelocity)
{
	// note: velocity u_{i+1/2} is practically stored at i+1
}

// Problem 3
void ExAdvectWithSemiLagrange(int xRes, int yRes, double dt,double* xVelocity, double* yVelocity, double *density, double *densityTemp,double *xVelocityTemp,double *yVelocityTemp)
{
	// note: velocity u_{i+1/2} is practically stored at i+1
}
