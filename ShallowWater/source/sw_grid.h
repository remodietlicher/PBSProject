/************************************************************
Sw_grid.h:
**********

2D quadrilateral grid:
* Sw_grid(int _xRes, int _yRes, float _dx)
    xRes : resolution in x direction
    yRes : resolution in y direction
    dx   : size of quadratic cell in x and y direction
************************************************************/

#ifndef _SWGRID_
    #define _SWGRID_
#include <cassert>

enum FIELDNAME { ETA, GROUND, HEIGHT, VELX, VELY };

class Sw_grid
{
public:
    Sw_grid(int _xRes, int _yRes, float _dx);
    void switchOldNew(FIELDNAME QUANTITY);

    const int xRes, yRes;
    const float dx;
    float **oldFields, **newFields;
private:

    void switch_ptr(float *&one, float *&two);

    float *eta, *ground, *height, *velx, *vely;
    float *eta_new, *ground_new, *height_new, *velx_new, *vely_new;
    float *normal;
//    int nFields;
};
#endif
