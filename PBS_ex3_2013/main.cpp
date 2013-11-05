// 2d fluid example, glut viewer

#include <iostream>
#include <string>
#include "util.h"
#include "fluid2d.h"

// toggle velocity display
// 0 = only show density
// 1 = show density, + velocity
// 2 = show density, + pressure
// 3 = show density, + divergence
int showFields = 0;

using namespace std;

// simulation variables
FLUID_2D* fluid = NULL;
bool pauseFlag = true;

// drawing variables
int width  = 800;
int height = 600;
double camera[] = { 0.0f, 1.0f, 0.0f, 1.0f };

// interaction
const double kernel[9] = {.5,.5,.5,.5,1,.5,.5,.5,.5};
const double sf = 4;
int oldx=-1, oldy=-1;

void reshapeGlut(int w, int h)
{
	if (h == 0) h = 1;

	glViewport(0, 0, w, h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(camera[0], camera[1], camera[2], camera[3], -10.0f, 10.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}


static void drawFields(const int xRes, const int yRes, const double* red, const double* green, const double* blue)
{
	glPushMatrix();
	glScalef(1.0f / (double)(xRes - 1), 1.0f / (double)(yRes - 1), 1.0f);

	double redMax = 0.0f;
	double greenMax = 0.0f;
	double blueMax = 0.0f;
	double redMin = ABS(red[1+xRes]);
	double greenMin = ABS(green[1+xRes]);
	double blueMin = ABS(blue[1+xRes]);
	for (int i = 0; i < xRes * yRes; i++)
	{
		redMax = ABS(red[i]) > redMax ? ABS(red[i]) : redMax;
		greenMax = ABS(green[i]) > greenMax ? ABS(green[i]) : greenMax;
		blueMax = ABS(blue[i]) > blueMax ? ABS(blue[i]) : blueMax;

		redMin = ABS(red[i]) < redMin ? ABS(red[i]) : redMin;
		greenMin = ABS(green[i]) < greenMin ? ABS(green[i]) : greenMin;
		blueMin = ABS(blue[i]) < blueMin ? ABS(blue[i]) : blueMin;
	}
	redMax = 1.0f / (redMax-redMin);
	greenMax = 1.0f / (greenMax-greenMin);
	blueMax = 1.0f / (blueMax-blueMin);

	for (int y = 0; y < yRes - 1; y++)
		for (int x = 0; x < xRes - 1; x++)
		{
			int index = x + y * xRes;

			glBegin(GL_TRIANGLE_FAN);
			double SWred = (ABS(red[index])-redMin) * redMax;
			double SEred = (ABS(red[index + 1])-redMin) * redMax;
			double NWred = (ABS(red[index + xRes])-redMin) * redMax;
			double NEred = (ABS(red[index + xRes + 1])-redMin) * redMax;
			double averagered = (SWred + SEred + NWred + NEred) * 0.25f;

			double SWgreen = (ABS(green[index])-greenMin) * greenMax;
			double SEgreen = (ABS(green[index + 1])-greenMin) * greenMax;
			double NWgreen = (ABS(green[index + xRes])-greenMin) * greenMax;
			double NEgreen = (ABS(green[index + xRes + 1])-greenMin) * greenMax;
			double averagegreen = (SWgreen + SEgreen + NWgreen + NEgreen) * 0.25f;

			double SWblue = (ABS(blue[index])-blueMin) * blueMax;
			double SEblue = (ABS(blue[index + 1])-blueMin) * blueMax;
			double NWblue = (ABS(blue[index + xRes])-blueMin) * blueMax;
			double NEblue = (ABS(blue[index + xRes + 1])-blueMin) * blueMax;
			double averageblue = (SWblue + SEblue + NWblue + NEblue) * 0.25f;

			glColor4f(averagered, averagegreen, averageblue, 1.0f);
			glVertex3f(x + 0.5f, y + 0.5f, 0.0f);
			glColor4f(SWred, SWgreen, SWblue, 1.0f);
			glVertex3f(x,     y, 0.0f);
			glColor4f(NWred, NWgreen, NWblue, 1.0f);
			glVertex3f(x    , y + 1, 0.0f);
			glColor4f(NEred, NEgreen, NEblue, 1.0f);
			glVertex3f(x + 1, y + 1, 0.0f);
			glColor4f(SEred, SEgreen, SEblue, 1.0f);
			glVertex3f(x + 1, y, 0.0f);
			glColor4f(SWred, SWgreen, SWblue, 1.0f);
			glVertex3f(x,     y, 0.0f);
			glEnd();
		}
	glPopMatrix();
}

// obtain idx into field from mouse position
bool unproject(int x, int y, int &rx, int &ry)
{
	double modelview[16], projection[16], ox,oy,oz;
    int viewport[4];
	float z;
    glGetDoublev( GL_PROJECTION_MATRIX, projection );
	glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
	glGetIntegerv( GL_VIEWPORT, viewport );
	glReadPixels( x, viewport[3]-y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z );
	gluUnProject( x, viewport[3]-y, z, modelview, projection, viewport, &ox,&oy,&oz );
	
	rx=(int)(ox*(fluid->xRes()-1)+.5);
	ry=(int)(oy*(fluid->yRes()-1)+.5);
	return !(rx<0 || rx>=fluid->xRes() || ry<0 || ry>=fluid->yRes());
}

void displayCallback()
{
	// setup for draw
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(camera[0], camera[1], camera[2], camera[3], -10.0f, 10.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glClear(GL_COLOR_BUFFER_BIT);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	switch(showFields) {
	case 0:
		drawFields(fluid->xRes(), fluid->yRes(), 
				fluid->density(), fluid->density(), fluid->density()); break;
	case 1:
		drawFields(fluid->xRes(), fluid->yRes(), 
				fluid->density(), fluid->xVelocity(), fluid->yVelocity()); break;
	case 2:
		drawFields(fluid->xRes(), fluid->yRes(), 
				fluid->density(), fluid->divergence(), fluid->divergence()); break;
	case 3:
		drawFields(fluid->xRes(), fluid->yRes(), 
				fluid->density(), fluid->divergence(), fluid->divergence()); break;
	}

	glutSwapBuffers();
}

void updateTitle()
{
	string title = "Exercise 4 - show density";
	switch (showFields)
	{
		case 0: break;
		case 1: title+="(r) + velocity(g,b)";break;
		case 2: title+="(r) + pressure(b)";break;
		case 3: title+="(r) + divergence(b)";break;
	}
	if (!pauseFlag) title+=" - PAUSED, press space";
	glutSetWindowTitle(title.c_str());
}

void keybCallback(unsigned char key, int x, int y)
{
	switch(key) 
	{
		// quit
		case 'q':
			exit(0);
			break;
		
		case 'd':
			showFields = (showFields+1) % 4;
			updateTitle();
			break;

		case ' ': 
			pauseFlag = !pauseFlag;
			updateTitle();
			break;
	}

	glutPostRedisplay();
}

void mouseCallback(int button, int state, int x, int y)
{
	unproject(x,y,oldx,oldy);
}

void motionCallback (int x,int y)
{
	int ox,oy;
	if (unproject(x,y,ox,oy) && oldx>=0)
	{
		for (int i=-1,ik=0;i<1;i++)
			for (int j=-1;j<1;j++,ik++)
			{
				if (ox+i<0 || ox+i>=fluid->xRes() || oy+j<0 || oy+j>=fluid->yRes()) continue;
				int idx=ox+i+fluid->xRes()*(oy+j);
				fluid->xForce()[idx] += sf * kernel[ik] * (ox-oldx);
				fluid->yForce()[idx] += sf *   kernel[ik] * (oy-oldy);				
			}	
	}
	printf("%d %d\n",ox-oldx,oy-oldy);
	oldx=ox;
	oldy=oy;
}

void idleFunc()
{ 
	if (pauseFlag)
	{
		static int frame = 0;

		// bottom inlet
		fluid->addDensity(.45,.55, .1,.15 );

		// centered for debugging
		//fluid->addDensity(.45,.55, .45,.55 );

		fluid->step(); 
		frame++;

		// force exit after a while
		if (frame == 10000) 
			exit(0);
	}

	glutPostRedisplay();
}


// setup glut 
int main(int argc, char **argv)
{
	// setup
	int res = SIM_RES;

	// print config
	cout << "2D Fluid test, res:"<<res<< endl << endl;
	cout << "Usage: start/stop='space', quit='q', toggle display='d'" << endl << endl;
	fluid = new FLUID_2D(res, res);

#if DO_PNG_OUT==1
	// unix-based systems only
	system("mkdir " OUTPATH);
	
	// delete any existing jpgs
	string path = string("./" OUTPATH);
	string remove = string("rm -f ") + path + string("/*.png");
	system(remove.c_str());
#endif

	// intialize GLUT window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(50, 50);
	glutInitWindowSize(width, height);
	glutCreateWindow("fluid2d");

	// init GLUT callbacks
	glutDisplayFunc(displayCallback);
	glutKeyboardFunc(keybCallback);
	glutMouseFunc(mouseCallback);
	glutMotionFunc(motionCallback);
	glutIdleFunc(idleFunc);  
	reshapeGlut(width, height); 

	// init GLUT drawing params
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	updateTitle();

	glutMainLoop();

	// clean up
	delete fluid;
	return 0;
}

