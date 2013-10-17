/*=====================================================================================*/
/*! \file		MeshViewer.h
	\author		peterkau
	\brief		Declaration of class MeshViewer
 */
/*=====================================================================================*/

#ifndef MESH_VIEWER_H
#define MESH_VIEWER_H

class FEMMesh;

//! Helper class for displaying an \c FEMMesh using OpenGL/GLUT
class MeshViewer
{
public:
	MeshViewer(int argc, char *argv[]);

	/*! Displays the solution \c solution for mesh \c mesh. */
	void VisualizeSolution(const FEMMesh &mesh, const vector<float> &solution);
	/*! Displays the error \c erro for mesh \c mesh. */
	void VisualizeError(const FEMMesh &mesh, const vector<float> &error);

private:
	//! GLUT display callback
	static void display();

	//! GLUT reshape callback
	static void reshape(int w, int h);

	//! Converts a HSV color to RGB
	static void HSV2RGB(float h, float s, float v, float &r, float &g, float &b);

	//! Creates geometry that represents the solution for mesh \c mesh
	static void createSolutionGeometry(const FEMMesh &mesh, const vector<float> &solution, float maxValue);
};

#endif
