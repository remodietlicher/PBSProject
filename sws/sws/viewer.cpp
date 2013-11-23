#include "viewer.h"
#include <stdlib.h>

SWViewer::SWViewer(SWSolver *sws){
	this->sws = sws;
	this->normals = new Vec3f[sws->getXRes()*sws->getYRes()];
	this->computedNormals = false;
}

void SWViewer::cleanup(){
	delete normals;
}

void SWViewer::computeNormals(){
	if(computedNormals) return;
	std::vector<float>& h = sws->getHeightMap();

//	Vec3f* normals2 = new Vec3f[sws->getXRes()*sws->getYRes()];
	for(int i=0; i<sws->getXRes(); i++)
		for(int j=0; j<sws->getYRes(); j++){
			Vec3f sum(0.0f, 0.0f, 0.0f);

			Vec3f out, in, left, right;
			if(i>0) {
				out = Vec3f(0.0f, h[V_INDEX(i-1, j)] - h[V_INDEX(i, j)], -1.0f);
			}
			if(i<sws->getXRes()-1){
				in = Vec3f(0.0f, h[V_INDEX(i+1, j)] - h[V_INDEX(i, j)], 1.0f);
			}
			if(j>0){
				left = Vec3f(-1.0f, h[V_INDEX(i, j-1)] - h[V_INDEX(i, j)], 0.0f);
			}
			if(i<sws->getYRes()-1){
				right = Vec3f(1.0f, h[V_INDEX(i, j+1)] - h[V_INDEX(i, j)], 0.0f);
			}

			if(j>0 && i>0){
				sum = sum + out.cross(left).normalize();
			}
			if(j>0 && i<sws->getXRes()-1){
				sum = sum + left.cross(in).normalize();
			}
			if(j<sws->getYRes()-1 && i<sws->getXRes()-1){
				sum = sum + in.cross(right).normalize();
			}
			if (j < sws->getYRes()-1 && i > 0) {
			sum = sum + right.cross(out).normalize();
			}

			normals[V_INDEX(i, j)] = sum;
		}
		computedNormals = true;
}

void SWViewer::handleKeypress(unsigned char key, int x, int y) {
	switch (key) {
		case 27: //Escape key
			exit(0);
	}
}

void SWViewer::initRendering() {
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_NORMALIZE);
	glShadeModel(GL_SMOOTH);
}

void SWViewer::handleResize(int w, int h) {
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (double)w / (double)h, 1.0, 200.0);
}

void SWViewer::drawScene() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -10.0f);
//	glRotatef(30.0f, 1.0f, 0.0f, 0.0f);
//	glRotatef(-_angle, 0.0f, 1.0f, 0.0f);
	
	GLfloat ambientColor[] = {0.4f, 0.4f, 0.4f, 1.0f};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientColor);
	
	GLfloat lightColor0[] = {0.6f, 0.6f, 0.6f, 1.0f};
	GLfloat lightPos0[] = {-0.5f, 0.8f, 0.1f, 0.0f};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor0);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
	
	glutSwapBuffers();
}

void SWViewer::createSolutionGeometry(){
	float scale = 5.0f / (sws->getXRes()-1);
	glScalef(scale, scale, scale);
	glTranslatef(-(float)(sws->getXRes()-1)/2,
				 0.0f,
				 -(float)(sws->getXRes()-1)/2);
	
	glColor3f(0.3f, 0.9f, 0.0f);
	std::vector<float>& h = sws->getHeightMap();
	for(int i = 0; i < sws->getXRes()-1; i++) {
		//Makes OpenGL draw a triangle at every three consecutive vertices
		glBegin(GL_TRIANGLE_STRIP);
		for(int j = 0; j < sws->getYRes(); j++) {
			Vec3f normal = getNormal(j, i);
			glNormal3f(normal(0), normal(1), normal(2));
			glVertex3f(j, h[V_INDEX(j, i)], i);
			normal = getNormal(j, i+1);
			glNormal3f(normal(0), normal(1), normal(2));
			glVertex3f(j, h[V_INDEX(j, i+1)], i+1);
		}
		glEnd();
	}
}

void SWViewer::viewScene(int argc, char** argv){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(400, 400);
	
	glutCreateWindow("Terrain - videotutorialsrock.com");
	initRendering();
	
	glutDisplayFunc(drawScene);
	glutKeyboardFunc(handleKeypress);
	glutReshapeFunc(handleResize);
//	glutTimerFunc(25, this->update, 0);
	createSolutionGeometry();
	
	glutMainLoop();
}

