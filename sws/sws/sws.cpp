#include "sws.h"
#include "vec3f.h"
#include <iostream>
#include <algorithm>
#include <fstream>


SWSolver::SWSolver(int xRes, int yRes, float xSize, float ySize, float dt) {
	// initialize defaults
	a_ext[0] = 0;
	a_ext[1] = 0;
	v_ext[0] = 0;
	v_ext[1] = 0;
	g = 9.81f;
	ground.resize(xRes*yRes);	// initializes all values to zero
	eta.resize(xRes*yRes);
	eta.assign(xRes*yRes, 1);	// put some default values
	height.resize(xRes*yRes);
	vel_x.resize(xRes*yRes);
//	vel_x.assign(xRes*yRes, 5.0f);
	vel_y.resize(xRes*yRes);
//	vel_x.assign(xRes*yRes, 5.0f);
	// rest
	this->res[0] = xRes;
	this->res[1] = yRes;
	this->xSize = xSize;
	this->ySize = ySize;
	this->dt = dt;
	this->dx[0] = xSize/xRes;
	this->dx[1] = ySize/yRes;

	calculateHeight();
}

float SWSolver::interpolate(std::vector<float> &array, float x, float y) {
	const int I = (int)x;
	const int J = (int)y;
	const float s1 = x - I;
	const float s0 = 1.f - s1;
	const float t1 = y - J;
	const float t0 = 1.f-t1;
	return  s0*(t0*array[INDEX(I, J)]+t1*array[INDEX(I, J+1)] )+s1*(t0*array[INDEX(I+1, J)] +t1*array[INDEX(I+1, J+1)]);
}

void SWSolver::advect(int QUANTITY){
	float v[2];

	std::vector<float> temp;
	std::vector<float> *array = nullptr;
	temp.resize(res[0]*res[1]);

	for(int i=1; i<res[0]-1; i++)
		for(int j=1; j<res[1]-1; j++){
			// set external velocities present at all times
			v[0] = v_ext[0];
			v[1] = v_ext[1];
			switch(QUANTITY){
			case 0: // eta (note that this is at the center of a cell)
				array = &eta;
				v[0] += (vel_x[INDEX(i, j)] + vel_x[INDEX(i+1, j)])*0.5f;
				v[1] += (vel_y[INDEX(i, j)] + vel_y[INDEX(i, j+1)])*0.5f;
				break;
			case 1: // vel_x (on the left)
				array = &vel_x;
				v[0] += vel_x[INDEX(i, j)];
				v[1] += (vel_y[INDEX(i, j)] + vel_y[INDEX(i-1, j)] + vel_y[INDEX(i-1, j+1)] + vel_y[INDEX(i, j+1)])*0.25f;
				break;
			case 2: // vel_y (on the bottom)
				array = &vel_y;
				v[0] += (vel_x[INDEX(i, j)] + vel_x[INDEX(i-1, j)] + vel_x[INDEX(i-1, j+1)] + vel_x[INDEX(i, j+1)])*0.25f;
				v[1] += vel_x[INDEX(i, j)];
				break;
			default: // error
				exit(1);
			}

			// backtrace position
			float srcpi = (float)i - v[0] * dt * 1/dx[0];
			float srcpj = (float)j - v[1] * dt * 1/dx[1];

			// clamp range of accesses (dont access boundaries (at 0 and res-1))
			if(srcpi<0.) srcpi = 1.f;
			if(srcpj<0.) srcpj = 1.f;
			if(srcpi>res[0]-1.) srcpi = res[0]-2.f;
			if(srcpj>res[1]-1.) srcpj = res[1]-2.f;

			// interpolate source value
			temp[INDEX(i, j)] = interpolate(*array, srcpi, srcpj);
		}
		for (int i = 0; i < res[0]; i++)
			for (int j = 0; j < res[1]; j++){
				(*array)[INDEX(i, j)] = temp[INDEX(i, j)];
	}
}

void SWSolver::updateHeight(){
	//std::cout << "update height:" << std::endl;
	for (int i = 1; i < res[0] - 1; i++) {
		//std::cout << std::endl << i << ": ";
		for (int j = 1; j < res[1] - 1; j++){

			float oldEta = eta[INDEX(i, j)];

			eta[INDEX(i, j)] += -eta[INDEX(i, j)] * dt * ((vel_x[INDEX(i + 1, j)] - vel_x[INDEX(i, j)]) / dx[0] + (vel_y[INDEX(i, j + 1)] - vel_y[INDEX(i, j)]) / dx[1]);
			if (eta[INDEX(i, j)] < ground[INDEX(i, j)]) eta[INDEX(i, j)] = ground[INDEX(i, j)];
			height[INDEX(i, j)] = eta[INDEX(i, j)] + ground[INDEX(i, j)];

			float newEta = eta[INDEX(i, j)];

			//std::cout << newEta - oldEta<<" ";
		}
		
	}
	//std::cout <<std::endl;
}

void SWSolver::updateVelocity(){
	for(int i=2; i<res[0]-1; i++)
		for(int j=1; j<res[1]-1; j++){
			vel_x[INDEX(i, j)] += (g*((height[INDEX(i-1, j)]-height[INDEX(i, j)])/dx[0])+a_ext[0])*dt;
		}
	for(int i=1; i<res[0]-1; i++)
		for(int j=2; j<res[1]-1; j++){
			vel_y[INDEX(i, j)] += (g*((height[INDEX(i, j-1)]-height[INDEX(i, j)])/dx[0])+a_ext[0])*dt;
		}
}

void SWSolver::setBoundary(){
	for(int i=0; i<res[0]; i++){
		height[INDEX(i, 0)] = height[INDEX(i, 1)];				// lower boundary
//		height[INDEX(i, 0)] = 7.0;				// lower boundary
		height[INDEX(i, res[0]-1)] = height[INDEX(i, res[0]-2)];	// upper boundary
	}
	for(int j=0; j<res[1]; j++){
		height[INDEX(0, j)] = height[INDEX(1, j)];				// left boundary
//		height[INDEX(0, j)] = 3.0;				// left boundary
		height[INDEX(res[1]-1, j)] = height[INDEX(res[1]-2, j)];	// right boundary
	}
}

void SWSolver::advanceTimestep(){
		//std::cout << "advecting eta..." << std::endl;
		advect(ETA);
		//std::cout << "advecting velocity_x..." << std::endl;
		advect(VELOCITY_X);
		//std::cout << "advecting velocity_y..." << std::endl;
		advect(VELOCITY_Y);

		//std::cout << "updating heights..." << std::endl;
		updateHeight();

		//std::cout << "updating velocities..." << std::endl;
		updateVelocity();

		//std::cout << "setting boundaries..." << std::endl;
		setBoundary();
}

void SWSolver::copyHeights(float* buffer)
{
	for (int i = 1; i < res[0] - 1; i++) {
		//std::cout << std::endl << i << ": ";
		for (int j = 1; j < res[1] - 1; j++){

			float h = height[INDEX(i, j)];
			buffer[INDEX(i, j)] = h;
		}
	}
}

const std::vector<float> SWSolver::getHeightMap(){
	return height;
}

void SWSolver::setEta(std::vector<float> eta){
	this->eta = eta;
	calculateHeight();
}
void SWSolver::setGround(std::vector<float> ground){
	this->ground = ground;
	calculateHeight();
}
void SWSolver::setVelocities(std::vector<float> vel_x, std::vector<float> vel_y){
	this->vel_x = vel_x;
	this->vel_y = vel_y;
}
void SWSolver::setExternalVelocities(float v_ext_x, float v_ext_y){
	this->v_ext[0] = v_ext_x;
	this->v_ext[1] = v_ext_y;
}
void SWSolver::setExternalAccelerations(float a_ext_x, float a_ext_y){
	this->a_ext[0] = a_ext_x;
	this->a_ext[1] = a_ext_y;
}
void SWSolver::calculateHeight(){
	for(int i=0; i<res[0]; i++)
		for(int j=0; j<res[1]; j++){
			height[INDEX(i, j)] = eta[INDEX(i, j)] + ground[INDEX(i, j)];
		}
}

float SWSolver::getXSize(){
	return xSize;
}
float SWSolver::getYSize(){
	return ySize;
}
int SWSolver::getXRes(){
	return res[0];
}
int SWSolver::getYRes(){
	return res[1];
}

SWRBSolver::SWRBSolver(int xRes, int yRes, float xSize, float ySize, float dt, Box *b) : 
	SWSolver(xRes, yRes, xSize, ySize, dt),
	box(b),
	alpha(1),
	rbs(b)
{
	displ_old.resize(xRes*yRes);
	displ_new.resize(xRes*yRes);
}

Box* SWRBSolver::getBody(){
	return this->box;
}

void SWRBSolver::advanceTimestep(){
	SWSolver::advanceTimestep();

	cout << "handling body interaction..." << endl;

	handleBodyInteraction();
}

void SWRBSolver::estimateIndices(Vector3f vertices[8], int &x_min, int &x_max, int &y_min, int &y_max){
//	std::vector<Vector3f> convHull;

	// for more accurate estimate of indices
/*	convHull = getConvexHull8XY(vertices);
	for(unsigned int i=0; i<convHull.size(); i++)
		cout << convHull[i][0] << " " << convHull[i][1] << " " << convHull[i][2] << endl;
*/

	// calculate bounding box on grid for a rough estimate of grid points
	// get points with x_max and x_min
	bubbleSortVert(0, vertices, 8);
	x_min = int(vertices[0][0]*res[0]);	// always round down for minimal index
	x_max = int(vertices[7][0]*res[0]+1);	// always round up for maximal index
	bubbleSortVert(1, vertices, 8);
	y_min = int(vertices[0][1]*res[1]);
	y_max = int(vertices[7][1]*res[1]+1);
}

void SWRBSolver::handleBodyInteraction(){
	// calculate vertices of the box (com + 1/2*(x0 + y0 + z0): +++, ++-, +--, ...)
	Vector3f vertices[8];
	int cnt = 0;
	for(int i=0; i<2; i++)
		for(int j=0; j<2; j++)
			for(int l=0; l<2; l++){
				vertices[cnt] = box->x + (box->x0*pow(-1, i)+box->y0*pow(-1, j)+box->z0*pow(-1, l))*0.5f;
				cnt++;
			}

	int x_min, x_max, y_min, y_max;
	estimateIndices(vertices, x_min, x_max, y_min, y_max);
	// note that here the iteration i depends on i-1 (when considering parallelization)
	// because the displacement depends on the height which is changed in every iteration
	Vector3f r;
	std::vector<Vector3f> positions;
	std::vector<Vector3f> forces;
	bool isIntersecting;
	for(int i=x_min; i<x_max; i++)
		for(int j=y_min; j<y_max; j++){
			isIntersecting = calculateDisplacement(i, j, displ_new[INDEX(i, j)], r);
			// handle body --> water
			// set height of neighbouring cells
			float dh = alpha*(displ_old[INDEX(i, j)]-displ_old[INDEX(i, j)])*0.25;
			height[INDEX(i+1, j)] += dh;
			height[INDEX(i, j+1)] += dh;
			height[INDEX(i-1, j)] += dh;
			height[INDEX(i, j-1)] += dh;

			// handle water --> body
			if(isIntersecting){
				float hb = 1;  // what is hb in Thue07 eq (11)? 
				float buoyancy = g*res[0]*res[1]*displ_new[INDEX(i, j)]*box->mass/(box->x0.length()+box->y0.length()+box->z0.length());
				forces.push_back(Vector3f(hb*vel_x[INDEX(i, j)], hb*vel_y[INDEX(i, j)], buoyancy));
				positions.push_back(r);
			}
		}

	// add gravity
	forces.push_back(Vector3f(0.0f, 0.0f, -9.81f));
	positions.push_back(box->x);

	rbs.advanceTimestep(dt, forces, positions);
}

//! b is the "height" of the displacement and r is the position on the surface of the body where the intersection happens
// r is later on used to calculate the torque excerted on the body
bool SWRBSolver::calculateDisplacement(int i, int j, float &displ, Vector3f &r){
	Vector3f P_line((float(i)+0.5f)/res[0], (float(j)+0.5f)/res[0], 0.0f);
	Vector3f V(0.0f, 0.0f, 1.0f);

	Vector3f vertex_bl(box->x - (box->x0 + box->y0 + box->z0)*0.5f);
	Vector3f vertex_tr(box->x + (box->x0 + box->y0 + box->z0)*0.5f);

	// define planes to make things easier
	Vector3f planeVector1[6];
	Vector3f planeVector2[6];
	Vector3f planePoint[6];
	planeVector1[0] = box->x0;
	planeVector2[0] = box->y0;
	planeVector1[1] = box->x0;
	planeVector2[1] = box->z0;
	planeVector1[2] = box->z0;
	planeVector2[2] = box->y0;
	for(int l=0; l<3; l++){
		planePoint[l] = vertex_bl;
		planeVector1[l+3] = -planeVector1[l];
		planeVector2[l+3] = -planeVector2[l];
		planePoint[l+3] = vertex_tr;
	}

	int counter = 0;
	bool isIntersecting;
	Vector3f P[2];
	// loop through planes and look for intersections
	for(int l=0; l<6; l++){
		Vector3f N = planeVector1[l].cross(planeVector2[l]);
		if((V|N) == 0){
			continue;
		}
		float t = (planePoint[l]-P_line)|N/(V|N);
		P[counter] = P_line + t*V;

		// if the intersection point is further away from the chosen vertex than p1+p2 its not on the box
		cout << (P[counter]-planePoint[l]).length() << " " <<  (planeVector1[l]+planeVector2[l]).length() << endl;
		if((P[counter]-planePoint[l]).length() < (planeVector1[l]+planeVector2[l]).length()) counter++;
	}
	// case of intersection (note that the box is convex)
	// for convenience set P[0] to be the lower end (lower z) of the body and P[1] to be the upper end (P[0]<P[1] always)
	if(counter == 2){
		bubbleSortVert(2, P, 2);
		float t0 = P[0][2];
		float t1 = P[1][2];
		r = P[0];
		assert(t0<t1);							// just double checking bubblesort
		float h = height[INDEX(i, j)];
		isIntersecting = true;					// remember where the intersection starts
		if(t1 < h) displ = t1 - t0;				// body under water
		else if(t0 > h){
			 displ = 0.0f;						// body above water (NO intersection water-body)
			 isIntersecting = false;
		}
		else displ = h - t0;					// only h-t0 is under water
	}
	// case of no intersection
	else if(counter == 0){
		displ = 0.0f;
		isIntersecting = false;
	}
	// something weird happened!
	else
		cout << "gugugugussuusuuu" << endl;

	return isIntersecting;
}

void SWRBSolver::bubbleSortVert(int coord, Vector3f *A, int n){
	bool swapped;
	do{
		swapped = false;
		for(int i=0; i<n-1; i++){
			if(A[i][coord] > A[i+1][coord]){
				Vector3f tmp = A[i];
				A[i] = A[i+1];
				A[i+1] = tmp;
				swapped = true;
			}
		}
		n--;
	} 
	while(swapped == true);
}

bool isLeft(Vector3f a, Vector3f b, Vector3f c){
    return ((b[0] - a[0])*(c[1] - a[1]) - (b[1] - a[1])*(c[0] - a[0])) > 0;
}

std::vector<Vector3f> SWRBSolver::getConvexHull8XY(Vector3f* vertices){
	bubbleSortVert(1, vertices, 8);
	if(vertices[0][0] == vertices[1][0]){
		bubbleSortVert(0, vertices, 8);
	}
	Vector3f pointOnHull = vertices[0];
	Vector3f endpoint;
	std::vector<Vector3f> P;
	int i=0;
	do{
		P.push_back(pointOnHull);
		endpoint = vertices[0];
		for(int j=0; j<8; j++){
			if((endpoint == pointOnHull) || isLeft(P[i], endpoint, vertices[j]))
				endpoint = vertices[j];
		}
		i++;
		pointOnHull = endpoint;
	} while(endpoint != P[0]);

	return P;
}

void SWRBSolver::testSorting(){
	ofstream vertices_data_x;
	ofstream vertices_data_y;
	ofstream vertices_data;
	vertices_data_x.open("testdata//verticesx.txt");
	vertices_data_y.open("testdata//verticesy.txt");
	vertices_data.open("testdata//vertices.txt");
	Vector3f vertices[8] = {Vector3f(10, 1, 1), Vector3f(1, 1, 0), Vector3f(0, 7, 0), Vector3f(0, 1, 2), Vector3f(1, 1, 0), Vector3f(0, 1, 2), Vector3f(0, 2, 4), Vector3f(1, 2, 0)};
	for(int i=0; i<8; i++){
		vertices_data << vertices[i][0] << endl;
		vertices_data << vertices[i][1] << endl;
		vertices_data << vertices[i][2] << endl;
	}
		
	int x_min, x_max, y_min, y_max;
	estimateIndices(vertices, x_min, x_max, y_min, y_max);
	cout << x_min << " " << x_max << " " << y_min << " " << y_max << endl;

}
