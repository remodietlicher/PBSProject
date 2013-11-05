#include "util.h"
#include "fluid2d.h"

#include <sstream>
#include <stdlib.h>
#if _WIN32
#else
#include <sys/time.h>
#endif



FLUID_2D::FLUID_2D(int xRes, int yRes) : _xRes(xRes), _yRes(yRes)
{
	// set simulation consts
	_dt = 0.1f;
	_impulseStrength = 5000.0f;

	// adapt timestep to acceleration
	_dt = TIMESTEP * sqrt((xRes+yRes)*0.5);
	_totalSteps = 0;

	// allocate solver
	_solver = new SolverGaussSeidel(solverIterations, _dt, solverAccuracy);

	// allocate arrays
	_totalCells   = _xRes * _yRes;
	_divergence   = new double[_totalCells];
	_pressure     = new double[_totalCells];
	_xVelocity    = new double[_totalCells];
	_yVelocity    = new double[_totalCells];
	_xVelocityTemp = new double[_totalCells];
	_yVelocityTemp = new double[_totalCells];
	_xForce       = new double[_totalCells];
	_yForce       = new double[_totalCells];
	_vorticity    = new double[_totalCells];
	_density      = new double[_totalCells];
	_densityTemp   = new double[_totalCells];

	cout << " Memory used by FLUID_2D: " << 14 * _totalCells * 4 / 1e6 << " MB, timestep="<<_dt<<"\n";

	// init arrays
	for (int x = 0; x < _totalCells; x++)
	{
		_divergence[x]   = 0.0f;
		_pressure[x]     = 0.0f;
		_xVelocity[x]    = 0.0f;
		_yVelocity[x]    = 0.0f;
		_xVelocityTemp[x] = 0.0f;
		_yVelocityTemp[x] = 0.0f;
		_xForce[x]       = 0.0f;
		_yForce[x]       = 0.0f;
		_vorticity[x]    = 0.0f;
		_density[x]      = 0.0f;
		_densityTemp[x]   = 0.0f;
	}
}

FLUID_2D::~FLUID_2D()
{
	if (_solver)      delete _solver;
	if (_density)     delete[] _density;
	if (_densityTemp)  delete[] _densityTemp;
	if (_divergence)  delete[] _divergence;
	if (_pressure)    delete[] _pressure;
	if (_xVelocity)   delete[] _xVelocity;
	if (_yVelocity)   delete[] _yVelocity;
	if (_xVelocityTemp) delete[] _xVelocityTemp;
	if (_yVelocityTemp) delete[] _yVelocityTemp;
	if (_xForce)       delete[] _xForce;
	if (_yForce)       delete[] _yForce;
	if (_vorticity)    delete[] _vorticity;
}
  
// step simulation once
void FLUID_2D::step()
{
	long int startt = getTime();
	// add in new forces
	addBuoyancy();
	addForce();

	// remove divergence
	solvePressure();

	// advect everything
	if(1) {
		advectValues();
	}

	// dump original densities
	if(1) {
		string prefix = string(OUTPATH "/density.");
		if(_totalSteps%JPG_OUT_INTERVAL==0) {
			dumpNumberedImage(_totalSteps, prefix, _density, _xRes, _yRes);
		} 
	}

	// reset forces
	for (int i = 0; i < _totalCells; i++)
		_xForce[i] = _yForce[i] = 0.0f;

	_totalSteps++;

	long int stopt = getTime();
	std::cout <<"Step "<<_totalSteps<<" took "<<timeToString(stopt-startt)<<" \n";
}

void FLUID_2D::addForce()
{
  for (int i = 0; i < _totalCells; i++) {
    _xVelocity[i] += _dt * _xForce[i];
    _yVelocity[i] += _dt * _yForce[i];
  }
}

void FLUID_2D::reset()
{
  for (int x = 0; x < _totalCells; x++) {
    _density[x]      = 0.0f;
    _densityTemp[x]   = 0.0f;
    _divergence[x]   = 0.0f;
    _pressure[x]     = 0.0f;
    _xVelocity[x]    = 0.0f;
    _yVelocity[x]    = 0.0f;
    _xVelocityTemp[x] = 0.0f;
    _yVelocityTemp[x] = 0.0f;
    _xForce[x]       = 0.0f;
    _yForce[x]       = 0.0f;
    _vorticity[x]    = 0.0f;
  }
  _totalSteps = 0;
}

void FLUID_2D::solvePressure()
{
	// copy out the boundaries 
	setNeumannX(_xVelocity);
	setNeumannY(_yVelocity);
	setZeroY(_xVelocity);
	setZeroX(_yVelocity);

	computeDivergence();

	// solve Poisson equation
	copyBorder(_pressure);
	_solver->solvePoisson(_xRes, _yRes, _pressure, _divergence);

	correctVelocities();
}


void FLUID_2D::addBuoyancy()
{
	int x,y,index;
	double scaling = 64.0f / (double)_xRes;

	// calculate vorticity
	index = 0;
	for (y = 0; y < _yRes; y++)
		for (x = 0; x < _xRes; x++, index++) {
			_yForce[index] += 0.1f * _density[index] * scaling;
		}
}

void FLUID_2D::addDensity(double xMin, double xMax, double yMin, double yMax)
{
	// add fluid square 
	for (int y = (int)(yMin * _yRes); y < (int)(yMax * _yRes); y++)
		for (int x = (int)(xMin * _xRes); x < (int)(xMax * _xRes); x++)
		{
			get(_density,x,y) = 1.0f;
		}
}

void FLUID_2D::computeDivergence()
{
	// calculate divergence
	const double dx = 1.0f / _xRes;
	const double idx = 1. / (dx);
	const double dx2 = dx*dx;
	const double idtx = 1. / (_dt * dx);

	double absDiv = 0.;
	for (int y = 1; y < _yRes - 1; y++) {
		for (int x = 1; x < _xRes - 1; x++)
			{
				const int index = y*_xRes+x;
				const double xComponent = (_xVelocity[index + 1] - _xVelocity[index]) * idtx;
				const double yComponent = (_yVelocity[index + _xRes] - _yVelocity[index]) * idtx;
				_divergence[index] = -(xComponent + yComponent)*dx2;
				//_divergence[index] = -(xComponent + yComponent);
			}
		}
	for (int i=0;i<_xRes*_yRes;i++)
	{
		//if (fabs(_divergence[i])>1e-6)
		//	printf("%d %e\n",i,_divergence[i]);
	}
}


// call exercises
extern void ExSolvePoisson(int xRes, int yRes, int _iterations, double _accuracy, double* field, double* b);
extern void ExCorrectVelocities(int _xRes, int _yRes, double _dt, const double* _pressure, double* _xVelocity, double* _yVelocity);	
extern void ExAdvectWithSemiLagrange(int xRes, int yRes, double dt,double* xVelocity, double* yVelocity, double *density, double *densityTemp,double *xVelocityTemp,double *yVelocityTemp);

void SolverGaussSeidel::solvePoisson(int xRes, int yRes, double *field, double *b)
{
	ExSolvePoisson(xRes,yRes,_iterations,_accuracy,field,b);
}
void FLUID_2D::correctVelocities()
{
	ExCorrectVelocities(_xRes, _yRes, _dt, _pressure, _xVelocity, _yVelocity);
}

void FLUID_2D::advectValues()
{
	ExAdvectWithSemiLagrange(_xRes,_yRes,_dt,_xVelocity,_yVelocity,_density,_densityTemp,_xVelocityTemp,_yVelocityTemp);
}

