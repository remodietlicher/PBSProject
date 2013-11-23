#ifndef VIEWER_H
#define VIEWER_H

#include <iostream>
#include <stdlib.h>
#include "util.h"
#include "sws.h"

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define V_INDEX(i, j) (i)*sws->getXRes() + (j)

class SWViewer {
private:
	Vec3f* normals;
	bool computedNormals;
	SWSolver* sws;
public:
	SWViewer(SWSolver *sws);
	static void handleKeypress(unsigned char key, int x, int y);
	static void initRendering();
	static void handleResize(int w, int h);
	static void drawScene();
	static void update(int value);
	void viewScene(int argc, char** argv);
	Vec3f getNormal(int i, int j);
	void createSolutionGeometry();
private:
	void computeNormals();
	void cleanup();
};

#endif