


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
#include "exporter.h"

using namespace std;

//#define USE_OPENMP

#ifdef USE_OPENMP
#include <omp.h>
// thread sync flag
int data_flag = 0;
#endif

// size of the image to render
unsigned int g_width;
unsigned int g_height;

// image buffer
float* g_buffer;

SWSolver *solver = nullptr;
SWViewer *viewer = nullptr;
Exporter * exporter = nullptr;


//foward declarations
//void run_solver();


void loadScene() {
	g_width  = 1000;
	g_height = 1000;

	int res[2] = { 100, 100 };
	
	solver = new SWSolver(res[0], res[1], 1, 1, 0.0006f);
	viewer = new SWViewer(solver);
	exporter = new Exporter();

	//load start values
	std::vector<std::vector<float>> heightmaps;
	std::vector<float> initEta;
	initEta.resize(res[0] * res[1]);
	float initialEta = 10.0f;
	for (int i = 0; i < res[0]; i++)
	for (int j = 0; j < res[1]; j++){
		if (i < 4 * res[0] / 6 && i>2 * res[0] / 6 && j<4 * res[1] / 6 && j>2 * res[1] / 6) initEta[INDEX(i, j)] = initialEta * 2.0;
		else initEta[INDEX(i, j)] = initialEta;
	}

	solver->setEta(initEta);
}
// export the rendered image
void exportImage() {

	float *image;

	if ((image = (float *)malloc(3 * g_width * g_height * sizeof(float))) == NULL) {
		fprintf(stderr, "Failed to allocate memory for image\n");
		return;
	}


	/* Copy the image into our buffer */
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadBuffer(GL_BACK_LEFT);
	glReadPixels(0, 0, g_width, g_height, GL_RGB, GL_FLOAT, image);

	
	if (!exporter->exportImage(image, g_width, g_height)) {
		std::cerr << "Failed to export image!" << std::endl;
	}
	else {
		std::cout << "Image file exported." << std::endl;
	}

	delete image;
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

	g_width = w;
	g_height = h;
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (double)w / (double)h, 1.0, 200.0);
}

#ifdef USE_OPENMP
void run_solver()
{
	int tid = omp_get_thread_num();

	while (true) {
		// wait until  data have been consumed/displayed (spinning, no mutex in openMP)
#pragma omp flush (data_flag)
		while (data_flag == 1){
#pragma omp flush (data_flag)
		}
		
		// compute next time step
		solver->advanceTimestep();

		// update viewer  (copy buffer, compute normals) 
		viewer->update();

		// notify renderer 
#pragma omp flush (data_flag)
		data_flag = 1;
#pragma omp flush (data_flag)

		//cout << "timestep computed (Thread " <<tid<<")"<< endl;
	}
}
#endif

void drawScene() {
	static int frmCnt = 0;


#ifdef USE_OPENMP
	// wait for data (spinning, no mutex in openMP)
#pragma omp flush (data_flag)
	while (data_flag == 0){
#pragma omp flush (data_flag)
	}
#else
	// compute next time step
	solver->advanceTimestep();
	// update viewer  (copy buffer, compute normals) 
	viewer->update();
#endif

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

	//FIXME: geometry: x, y (grid) and z must match
	float scaleHeight = 1.0;
	float xRes = viewer->getXSize();
	float yRes = viewer->getYSize();

	glColor3f(0.1f, 0.1f, 0.8f);
	for (int z = 0; z < length - 1; z++) {
		//Makes OpenGL draw a triangle at every three consecutive vertices
		glBegin(GL_TRIANGLE_STRIP);
		for (int x = 0; x < width; x++) {
			Vec3f normal = viewer->getNormal(x, z);
			float height = scaleHeight * viewer->getHeight(x, z);
			glNormal3f(normal[0],normal[1], normal[2]);
			glVertex3f(x * xRes, height, z* yRes);

			normal = viewer->getNormal(x, z + 1);
			height = scaleHeight * viewer->getHeight(x, z+1);
			glNormal3f(normal[0], normal[1],normal[2]);
			glVertex3f(x* xRes, height, (z + 1) * yRes);
		}
		glEnd();
	}

	glutSwapBuffers();

	//use with debuger or or file system is getting flooded ..
	//exportImage();

#ifdef USE_OPENMP
	// notify solver that we are done
#pragma omp flush (data_flag)
	data_flag = 0;
#pragma omp flush (data_flag)

	int tid = omp_get_thread_num();
	//cout << "data consumed (Thread " << tid << ")" << endl;
#endif


	frmCnt++;
	if (frmCnt % 100 == 0) {
		cout << "number of rendered frames: " << frmCnt << endl;
	}
}

void update(int value) {
	
	glutPostRedisplay();
	glutTimerFunc(1, update, 0);
}

#ifdef USE_OPENMP
void start() 
{
#pragma omp parallel sections num_threads(2)
	{
#pragma omp section 
		{
			glutMainLoop();
		}

#pragma omp section 
		{
			run_solver();
		}
	}
}
#endif

int main(int argc, char** argv) {

	loadScene();

/*
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(g_width, g_height);

	glutCreateWindow("PBS HS2013 - shallow water");
	initRendering();


	glutDisplayFunc(drawScene);
	glutKeyboardFunc(handleKeypress);
	glutReshapeFunc(handleResize);
	glutTimerFunc(100, update, 0);

#ifdef USE_OPENMP
	start();
#else
	glutMainLoop();
#endif
*/
	SWRBSolver swrb(0, 0, 0.0f, 0.0f, 0.0f, Box(Vector3f(0, 0, 0), 0.0f, Vector3f(0, 0, 0), Vector3f(0, 0, 0), Vector3f(0, 0, 0)));

	swrb.testSorting();


	return 0;
}



