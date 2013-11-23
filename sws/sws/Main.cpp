


#include <iostream>
#include <stdlib.h>
#include <algorithm>

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include "GL/glut.h"
#endif

#include "sws.h"
#include "viewer.h"

using namespace std;

SWSolver *solver = nullptr;
SWViewer *viewer = nullptr;


void loadScene() {
	int res[2] = { 100, 100 };
	
	solver = new SWSolver(res[0], res[1], 1, 1, 0.001f);
	viewer = new SWViewer(solver);

	//load start values
	std::vector<std::vector<float>> heightmaps;
	std::vector<float> initEta;
	initEta.resize(res[0] * res[1]);
	for (int i = 0; i < res[0]; i++)
	for (int j = 0; j < res[1]; j++){
		if (i<4 * res[0] / 6 && i>2 * res[0] / 6 && j<4 * res[1] / 6 && j>2 * res[1] / 6) initEta[INDEX(i, j)] = 1.1f;
		else initEta[INDEX(i, j)] = 1.0f;
	}

	solver->setEta(initEta);
}

void cleanup() {
	//FXIME: clean up solver viewer
}

void handleKeypress(unsigned char key, int x, int y) {
	switch (key) {
	case 27: //Escape key
		cleanup();
		exit(0);
	}
}

void initRendering() {
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_NORMALIZE);
	glShadeModel(GL_SMOOTH);
}

void handleResize(int w, int h) {
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (double)w / (double)h, 1.0, 200.0);
}

void drawScene() {

	// compute next time step
	solver->advanceTimestep();

	// update normales 
	viewer->computeNormals();

	int width  = viewer->width();
	int length = viewer->length();

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -10.0f);
	glRotatef(30.0f, 1.0f, 0.0f, 0.0f);
	//glRotatef(-_angle, 0.0f, 1.0f, 0.0f);

	GLfloat ambientColor[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientColor);

	GLfloat lightColor0[] = { 0.6f, 0.6f, 0.6f, 1.0f };
	GLfloat lightPos0[] = { -0.5f, 0.8f, 0.1f, 0.0f };
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor0);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);

	float scale = 5.0f / max(width - 1, length - 1);
	glScalef(scale, scale, scale);
	glTranslatef(-(float)(width - 1) / 2,
		0.0f,
		-(float)(length - 1) / 2);

	glColor3f(0.3f, 0.9f, 0.0f);
	for (int z = 0; z < length - 1; z++) {
		//Makes OpenGL draw a triangle at every three consecutive vertices
		glBegin(GL_TRIANGLE_STRIP);
		for (int x = 0; x < width; x++) {
			Vec3f normal = viewer->getNormal(x, z);
			glNormal3f(normal[0], normal[1], normal[2]);
			glVertex3f(x,viewer->getHeight(x, z), z);

			normal = viewer->getNormal(x, z + 1);
			glNormal3f(normal[0], normal[1], normal[2]);
			glVertex3f(x, viewer->getHeight(x, z + 1), z + 1);
		}
		glEnd();
	}

	glutSwapBuffers();
}

void update(int value) {
	
	glutPostRedisplay();
	glutTimerFunc(25, update, 0);
}

int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(400, 400);

	glutCreateWindow("PBS HS2013 - shallow water");
	initRendering();

	loadScene();

	glutDisplayFunc(drawScene);
	glutKeyboardFunc(handleKeypress);
	glutReshapeFunc(handleResize);
	glutTimerFunc(25, update, 0);

	glutMainLoop();
	return 0;
}



