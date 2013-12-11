/************************************************************
Sw_grid.h:
**********

2D quadrilateral grid:
* Sw_grid(int _xRes, int _yRes, val _dx)
    xRes : resolution in x direction
    yRes : resolution in y direction
    dx   : size of quadratic cell in x and y direction

* val get(float field, int ix, int iy)
    use to READ values, ix/iy out of range -> boundary used

* val set(float field, int ix, int iy)
    use to WRITE values, ix/iy out of range -> error!

Example Usage:
    // set values
    for (int ix = 0; ix < grid.xRes; ++ix)
    for (int iy = 0; iy < grid.yRes; ++iy)
    {
        val value = (double)ix + 100.*(double)iy;
        grid.set(ETA, ix, iy) = some_value;
        // or:
        val *to_somewhere = &grid.set(ETA, ix, iy);
        *to_somewhere = some_value;
    }

    // get values
    for (int ix = -1; ix < grid.xRes+1; ++ix)
    for (int iy = -1; iy < grid.yRes+1; ++iy)
        printf("%02d/%02d : %f\n", ix, iy, grid.get(ETA, ix, iy));
    
************************************************************/

#include <cassert>
#include "sw_grid.h"

/* GRID */
Sw_grid::Sw_grid(int _xRes, int _yRes, float _dx) : xRes(_xRes), yRes(_yRes), dx(_dx) {
    
    int nFields = 5;

    int xySize = xRes*yRes;
    eta        = new float[xySize];
    ground     = new float[xySize];
    height     = new float[xySize];
    velx       = new float[xySize];
    vely       = new float[xySize];
    eta_new    = new float[xySize];
    ground_new = new float[xySize];
    height_new = new float[xySize];
    velx_new   = new float[xySize];
    vely_new   = new float[xySize];

    oldFields = new float*[nFields];

    oldFields[ETA]    = eta;
    oldFields[GROUND] = ground;
    oldFields[HEIGHT] = height;
    oldFields[VELX]   = velx;
    oldFields[VELY]   = vely;

    newFields = new float*[nFields];

    newFields[ETA]    = eta_new;
    newFields[GROUND] = ground_new;
    newFields[HEIGHT] = height_new;
    newFields[VELX]   = velx_new;
    newFields[VELY]   = vely_new;

    normal = new float[xySize];
}

void Sw_grid::switchOldNew(FIELDNAME QUANTITY){
    // for (int i = 0; i < nFields; ++i)
        switch_ptr(newFields[QUANTITY], oldFields[QUANTITY]);
}

void Sw_grid::switch_ptr(float *&one, float *&two){
    float *tmp;
    tmp = one;
    one = two;
    two = tmp;
}
