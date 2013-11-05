#ifndef FLUID_2D_H
#define FLUID_2D_H

#include <cmath>
#include <vector>
#include "util.h"

#include "solver.h"
#include "gauss_seidel.h"
#ifdef _WIN32
  #include "GLUT/glut.h"
#else
#ifdef APPLE
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif
#endif  

#include <iostream>

using namespace std;

// resolution of the simulation
const int SIM_RES = 23;

// global control defines
#define OUTPATH "./output/"
const int JPG_OUT_INTERVAL = 2;
const double TIMESTEP = 0.005;

// accuracy for gauss seidel solver
const double solverAccuracy = 1e-5;
// max. number of iterations
const int solverIterations = 1000;


// solver class
class FLUID_2D  
{
	public:
		FLUID_2D(int xRes, int yRes);

		// destructor
		virtual ~FLUID_2D();

		// accessors 
		double* density()   { return _density; }; 
		double* xVelocity() { return _xVelocity; };
		double* yVelocity() { return _yVelocity; };
		double* xForce() { return _xForce; };
		double* yForce() { return _yForce; };
		double* pressure()  { return _pressure; };
		double* divergence()  { return _divergence; };
		int xRes()         { return _xRes; };
		int yRes()         { return _yRes; };
		double dt()         { return _dt; };


		// Advance the simulation
		void step();

		// reset everything
		void reset();

		// add a density region 
		void addDensity(double xMin, double xMax, double yMin, double yMax);

	protected:
		int _xRes;              // x resolution of simulation
		int _yRes;              // y resolution of simulation
		int _totalCells;        // total cells in the simulation = _xRes * _yRes
		double _dt;              // simulation timestep
		double _impulseStrength; // strength of addImpulse(), used for debugging purposes
		int _totalSteps;        // total timesteps taken so far

		double* _density;        // current density field
		double* _densityTemp;     // previous density field
		double* _pressure;       // current pressure field - x in the Poisson solver
		double* _xVelocity;      // current x velocity
		double* _yVelocity;      // current y velocity
		double* _xVelocityTemp;   // previous x velocity
		double* _yVelocityTemp;   // previous y velocity
		double* _xForce;         // current x force
		double* _yForce;         // current y force
		double* _divergence;     // velocity divergence - b in the Poisson solver
		double* _vorticity;      // vorticity of the velocity field
		SolverBase* _solver;    // whatever Poisson solver you chose

		// solve the Poisson equation for the pressure correction
		void solvePressure();

		// add the x and y forces into the velocity field
		virtual void addForce();

		// add buoyancy into the velocity
		void addBuoyancy();


		// helper functions

		// compute divergence
		void computeDivergence();
		// correct velocities
		void correctVelocities();
		// advect values
		void advectValues();


		// set the x direction to periodic 
		void setPeriodicX(double* field);

		// set neumann boundaries
		void setNeumannX(double* field);
		void setNeumannY(double* field);
		void setZeroX(double* field);
		void setZeroY(double* field);
		void setZeroBorder(double* field) {
			setZeroX(field);
			setZeroY(field);
		};

		// copy out boundaries
		void copyBorderX(double* field);
		void copyBorderY(double* field);
		void copyBorder(double* field) {
			copyBorderX(field);
			copyBorderY(field);
		};

		// clamp a field to a min and max
		void clamp(double* field, double min, double max) {
			for (int x = 0; x < _xRes * _yRes; x++)
				field[x] = (field[x] > max) ? max : (field[x] < min) ? min : field[x];
		};

		// swap the old and new velocity fields
		void swapVelocity() {
			double* temp   = _xVelocity;
			_xVelocity    = _xVelocityTemp;
			_xVelocityTemp = temp;
			temp          = _yVelocity;
			_yVelocity    = _yVelocityTemp;
			_yVelocityTemp = temp;
		};

		// swap the old and new density fields
		void swapDensity() {
			double* temp = _density;
			_density = _densityTemp;
			_densityTemp = temp;
		};

		// access field in the solver
		// add index checks for debug mode
		template <class T > T& get(T* field, int i,int j) {
			return field[i + j * _xRes];
		}
};

#endif
