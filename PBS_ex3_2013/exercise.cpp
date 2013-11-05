#include "gauss_seidel.h"
#include "fluid2d.h"

#define GETF(x,y) field[(x) + (y) * xRes]
#define GETB(x,y) b[(x) + (y) * xRes]

// Problem 1
void ExSolvePoisson(int xRes, int yRes, int _iterations, double _accuracy, double* field, double* b)
{

	// double dt = fluid.dt;

	double *p_old = new double[xRes*yRes];
	double *p_new = new double[xRes*yRes];

	for (int y = 0; y < yRes; y++)
	for (int x = 0; x < xRes; x++){
		p_old[x + y*xRes] = 0.0;
		p_new[x + y*xRes] = GETF(x, y);
		// printf("%f\n", p_new[x + y*xRes]);
	}

	for (int iter = 0; iter < _iterations; ++iter)
	{
		double *p_old_ = p_old;
		p_old = p_new;
		p_new = p_old_;

		double residual = 0;

		// note that the boundaries are handles by the framework, so you iterations should be similar to:
		for (int y = 1; y < yRes - 1; y++) {
			for (int x = 1; x < xRes - 1; x++)
			{
				double d = GETB(x, y);
				double p_old_FD = p_old[x+1 + y*xRes] + p_old[x + (y+1)*xRes] + p_old[x-1 + y*xRes] + p_old[x + (y-1)*xRes];
				p_new[x + y*xRes] = 0.25 * ( d  + p_old_FD );

				double Ap_new = 4.0 * p_new[x + y*xRes] - (p_new[x+1 + y*xRes] + p_new[x + (y+1)*xRes] + p_new[x-1 + y*xRes] + p_new[x + (y-1)*xRes]);
				residual += fabs(d - Ap_new);
				// printf("b = %f\n", p_old[x+y*xRes]);
			}
		}

		residual /= (double)(xRes*yRes);


		//for your debugging, and ours, please add these prints after every iteration
		if(iter==_iterations-1)
			printf("Pressure solver: iter=%d , res=%f \n",iter, residual);
		if(residual<_accuracy) {
			printf("Pressure solver: iter=%d , converged \n",iter,residual);
			break; // optional
		}
	}

	for (int y = 1; y < yRes - 1; y++)
	for (int x = 1; x < xRes - 1; x++){
		GETF(x, y) = p_new[x + y*xRes];
	}


}

// Probelm 2
void ExCorrectVelocities(int _xRes, int _yRes, double _dt, const double* _pressure, double* _xVelocity, double* _yVelocity)
{
	double dx = 1.0/(double)_xRes;
	double dy = 1.0/(double)_yRes;
	// note: velocity u_{i+1/2} is practically stored at i+1
	for (int y = 0; y < _yRes; y++)
	for (int x = 0; x < _xRes; x++){
		_xVelocity[x+1 + y*_xRes]   -= _dt / dx * (_pressure[x+1 + y*_xRes] - _pressure[x + y*_xRes]);
		_yVelocity[x + (y+1)*_xRes] -= _dt / dy * (_pressure[x + (y+1)*_xRes] - _pressure[x + y*_xRes]);
	}
}

// Problem 3
void ExAdvectWithSemiLagrange(int xRes, int yRes, double dt, double* xVelocity, double* yVelocity, double *density, double *densityTemp, double *xVelocityTemp, double *yVelocityTemp)
{
	// note: velocity u_{i+1/2} is practically stored at i+1
	double dx[2];
	dx[0] = 1.0/(double)xRes;
	dx[1] = 1.0/(double)yRes;

	for (int y = 1; y < yRes-1; y++)
	for (int x = 1; x < xRes-1; x++){

		double adv[2];

		adv[0] = 0.5 * (xVelocity[x + y*xRes] + xVelocity[x+1 + y*xRes]);
		adv[1] = 0.5 * (yVelocity[x + y*xRes] + yVelocity[x + (y+1)*xRes]);

		double xg[] = {(double)x*dx[0], (double)y*dx[1]};

		double xp[2], weight[2]; // weight[2]: interpolation weight [0 to 1, where 0 takes x0 and 1 takes xp]
		int i0[2]; // index of bottom left gridpoint of the cell, where we interpolate
		for (int i = 0; i < 2; ++i)
		{
			xp[i]  = xg[i] - dt *adv[i];
			i0[i] = (int)( floor(xp[i]/dx[i]) );
			double x0 = dx[i] * (double)(i0[i]);
			weight[i] = ( xp[i] - x0 ) / dx[i];
			if(weight[i] < -0.000001 or weight[i] > 1.0000001)
			{
				printf("weight = %f\n", weight[i]);
				printf("dx = %f\n", dx[i]);
				printf("xp = %f\n", xp[i]);
				printf("x0 = %f\n", x0);
				exit(0);
			}
		}



		xVelocityTemp[x + y*xRes] = 0.0;
		yVelocityTemp[x + y*xRes] = 0.0;
		densityTemp[x + y*xRes]   = 0.0;

		int x0y0 = i0[0]   + ( i0[1]  )*xRes;
		int x1y0 = i0[0]+1 + ( i0[1]  )*xRes;
		int x0y1 = i0[0]   + ( i0[1]+1)*xRes;
		int x1y1 = i0[0]+1 + ( i0[1]+1)*xRes;
		xVelocityTemp[x + y*xRes] = (1.0-weight[1])*( (1.0-weight[0])*xVelocity[x0y0] + weight[0]*xVelocity[x1y0] ) + weight[1]*( (1.0-weight[0])*xVelocity[x0y1] + weight[0]*xVelocity[x1y1] );
		yVelocityTemp[x + y*xRes] = (1.0-weight[1])*( (1.0-weight[0])*yVelocity[x0y0] + weight[0]*yVelocity[x1y0] ) + weight[1]*( (1.0-weight[0])*yVelocity[x0y1] + weight[0]*yVelocity[x1y1] );
		densityTemp[x + y*xRes]   = (1.0-weight[1])*( (1.0-weight[0])*density[x0y0]   + weight[0]*density[x1y0] )   + weight[1]*( (1.0-weight[0])*density[x0y1]   + weight[0]*density[x1y1] );
	}
	for (int y = 0; y < yRes; y++)
	for (int x = 0; x < xRes; x++){
		xVelocity[x + y*xRes] = xVelocityTemp[x + y*xRes];
		yVelocity[x + y*xRes] = yVelocityTemp[x + y*xRes];
		density[x + y*xRes]   = densityTemp[x + y*xRes];
	}
}
