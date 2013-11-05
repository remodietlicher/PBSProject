#include "fluid2d.h"



void FLUID_2D::setPeriodicX(double* field)
{
	for (int y = 0; y < _yRes; y++)
	{
		int index = y * _xRes;
		field[index] = field[index + _xRes - 2];
		field[index + _xRes - 1] = field[index + 1];
	}
}

void FLUID_2D::copyBorderX(double* field)
{
	const int sx= _xRes;
	const int sy= _yRes;
	for (int y = 0; y < sy; y++)
	{
		int index = y * sx;
		field[index]             = field[index + 1];
		field[index + sx - 1] = field[index + (sx - 1) - 1];
	}
}

void FLUID_2D::copyBorderY(double* field)
{
	const int sx= _xRes;
	const int sy= _yRes;
	for (int x = 0; x < sx; x++)
	{
		int index = x;
		field[index] = field[index + sx];

		index = x + (sy - 1) * sx;
		field[index] = field[index - sx];
	}
}

void FLUID_2D::setNeumannX(double* field)
{
	for (int y = 0; y < _yRes; y++)
	{
		int index = y * _xRes;
		field[index]             = field[index + 2];
		field[index + _xRes - 1] = field[index + (_xRes - 1) - 2];
	}
}

void FLUID_2D::setNeumannY(double* field)
{
	for (int x = 0; x < _xRes; x++)
	{
		int index = x;
		field[index] = field[index + 2 * _xRes];

		index = x + (_yRes - 1) * _xRes;
		field[index] = field[index - 2 * _xRes];
	}
}


void FLUID_2D::setZeroX(double* field)
{
	const int sx= _xRes;
	const int sy= _yRes;
	for (int y = 0; y < sy; y++)
	{
		int index = y * sx;
		field[index]             = 0.0f;
		field[index + sx - 1] = 0.0f;
	}
}

void FLUID_2D::setZeroY(double* field)
{
	const int sx= _xRes;
	const int sy= _yRes;
	for (int x = 0; x < sx; x++)
	{
		int index = x;
		field[index] = 0.0f;

		index = x + (sy - 1) * sx;
		field[index] = 0.0f;
	}
}

