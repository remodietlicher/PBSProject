#ifndef VIEWER_H
#define VIEWER_H

#include <iostream>
#include <stdlib.h>
#include "vec3f.h"
#include "sws.h"

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include "GL/glut.h"
#endif

#define V_INDEX(i, j) (i)*sws->getXRes() + (j)

class SWViewer {
private:
	Vec3f* normals;
	float* heights;
	SWRBSolver* sws;
public:
	SWViewer(SWRBSolver *sws);
	int width();
	int length();
	float getXSize();
	float getYSize();

	void update();
	Vec3f getNormal(int i, int j);
	float getHeight(int i, int j);

private:
	void computeNormals();
	void cleanup();
};

#endif