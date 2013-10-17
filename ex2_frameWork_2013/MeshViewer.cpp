/*=====================================================================================*/
/*! \file		MeshViewer.cpp
	\author		peterkau
	\brief		Implementation of class MeshViewer
 */
/*=====================================================================================*/

#include <stdlib.h>
#include "SimpleFEMDefs.h"
#include "MeshViewer.h"
#include "FEMMesh.h"

#if defined(__APPLE__) || defined(__APPLE_CC__)
	#include <GLUT/glut.h>
#else
	#define GLUT_DISABLE_ATEXIT_HACK
	#include "GL/glut.h"	// Also includes GL/gl.h and GL/glu.h
#endif

MeshViewer::MeshViewer(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
}

void handle_keyboard(unsigned char key, int x, int y) {
	switch (key) {
		case 'q': 
		case 27: // ESC
			exit(0);
			break;
	}
}

void MeshViewer::VisualizeSolution(const FEMMesh &mesh, const vector<float> &solution)
{
	float maxValue = 0;
	for(int i = 0; i<solution.size(); i++)
		maxValue = std::max(solution[i],maxValue);
	
	glutCreateWindow("SimpleFEM");

	glutReshapeWindow(600, 600);

	glutDisplayFunc(display);

	glutReshapeFunc(reshape);

	glNewList(1, GL_COMPILE); 

	createSolutionGeometry(mesh, solution, maxValue);

	glEndList();

	glutKeyboardFunc(handle_keyboard);

	glutMainLoop();
}

void MeshViewer::VisualizeError(const FEMMesh &mesh, const vector<float> &error)
{
	int n = error.size();
	std::vector<float> abserror(error);
	for(int i=0; i<n; i++)
		abserror[i] = fabs(abserror[i]);

	VisualizeSolution(mesh, abserror);
}

void MeshViewer::display()
{
	// Clear the window
	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	// Map coordinates from [0,0] - [1,1] to [-1,-1] - [1,1]:
	glTranslatef(-1.0f, -1.0f, 0.0f);
	glScalef(2.0f, 2.0f, 1.0f);

	// Render the display list:
	glCallList(1);      

	glPopMatrix();

	// Bring the background buffer to the front:
	glutSwapBuffers();
}

void MeshViewer::reshape(int w, int h)
{
	glViewport(0, 0, w, h);
}

void MeshViewer::HSV2RGB(float h, float s, float v, float &r, float &g, float &b)
{
	// Get hue in [0,1):
	h /= 360.0;
	if(h == 1)
		h = 0;

	h *= 6.0;

	// Get used segment i, i in [0,5]:
	int i = floor(h);

	// Get the fractional part of h:
	float f = h-i;
	float p = v*(1-s);
	float q = v*(1-(s*f));
	float t = v*(1-(s*(1-f)));

	// Cases for every segment:
	switch (i)
	{
	case 0: r=v; g=t; b=p; break;
	case 1: r=q; g=v; b=p; break;
	case 2: r=p; g=v; b=t; break;
	case 3: r=p; g=q; b=v; break;
	case 4: r=t; g=p; b=v; break;
	case 5: r=v; g=p; b=q; break;
	}
}

void MeshViewer::createSolutionGeometry(const FEMMesh &mesh, const vector<float> &solution, float maxValue)
{

	// draw elements and solution
	glBegin(GL_TRIANGLES);

	// loop over elements
	for(size_t elID=0; elID<mesh.GetNumElements(); elID++)
	{   
		// get one
		const FEMElementTri &element = mesh.GetElement(elID);

		assert(element.GetNumElementNodes() == 3);

		// loop over nodes of element
		for(size_t elsize_t=0; elsize_t<3; elsize_t++)
		{
			// get global id
			size_t nodeID = element.GetGlobalNodeForElementNode(elsize_t);

			// get global position
			const Vector2 &pos = mesh.GetNodePosition(nodeID);
			float val = solution[nodeID] / maxValue;

			// convert value to color
			float s = 1.0f;
			float v = 1.0f;

			// Out of range -> black / red
			if(val < 0.0f ) {
				val = 0.0f;
				v = 0.0f;
			}
			if(val > 1.0f) {
				val = 1.0f;
				v = 359.0f;
			}

			float h = (1.0f - val) * 240.0f;

			float r, g, b;
			HSV2RGB(h, s, v, r, g, b);

			// draw
			glColor3f(r, g, b);
			glVertex3f(pos[0], pos[1], 0);
		}
	}

	glEnd();


	// overdraw black lines
	glColor3f(0., 0., 0.);
	glBegin(GL_LINES);

	for(size_t elID=0; elID<mesh.GetNumElements(); elID++)
	{   
		const FEMElementTri &element = mesh.GetElement(elID);

		for(size_t elsize_t=0; elsize_t<3; elsize_t++)
		{
			size_t nodeID1 = element.GetGlobalNodeForElementNode(elsize_t);
			size_t nodeID2 = element.GetGlobalNodeForElementNode((elsize_t+1)%3);

			const Vector2 &pos1 = mesh.GetNodePosition(nodeID1);
			const Vector2 &pos2 = mesh.GetNodePosition(nodeID2);

			glVertex3f(pos1[0], pos1[1], 0);
			glVertex3f(pos2[0], pos2[1], 0);
		}
	}

	glEnd();
}


