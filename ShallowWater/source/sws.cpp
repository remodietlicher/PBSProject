#include "sws.h"
#include "vec3f.h"
#include <iostream>
#include <math.h>


SWSolver::SWSolver() {
	// initialize defaults
	a_ext[0] = 0;
	a_ext[1] = 0;
	v_ext[0] = 0;
	v_ext[1] = 0;
	g = 9.81f;
	time = 0.0;

	// calculateHeight();
}

void SWSolver::setGrid(Sw_grid *_grid){
	grid = _grid;

	calculateHeight();
}

float SWSolver::interpolate(const float *array, float x, float y) {
	const int I = (int)x;
	const int J = (int)y;
	const float s1 = x - I;
	const float s0 = 1.f - s1;
	const float t1 = y - J;
	const float t0 = 1.f-t1;
	return  s0*(t0*array[INDEX(I, J)]+t1*array[INDEX(I, J+1)] )+s1*(t0*array[INDEX(I+1, J)] +t1*array[INDEX(I+1, J+1)]);
}

void SWSolver::advect(FIELDNAME QUANTITY){
	float v[2];

	const float *array = grid->oldFields[QUANTITY];
	const float *vel_x = grid->oldFields[VELX];
	const float *vel_y = grid->oldFields[VELY];

	float *noi  = grid->newFields[QUANTITY];

	float dx = grid->dx;
	int xRes = grid->xRes;
	int yRes = grid->yRes;

	for(int i=1; i<xRes-1; i++)
		for(int j=1; j<yRes-1; j++){
			// set external velocities present at all times
			v[0] = v_ext[0];
			v[1] = v_ext[1];
			switch(QUANTITY){
			case ETA: // eta (note that this is at the center of a cell)
				v[0] += (vel_x[INDEX(i, j)] + vel_x[INDEX(i+1, j)])*0.5f;
				v[1] += (vel_y[INDEX(i, j)] + vel_y[INDEX(i, j+1)])*0.5f;
				break;
			case VELX: // vel_x (on the left)
				v[0] += vel_x[INDEX(i, j)];
				v[1] += (vel_y[INDEX(i, j)] + vel_y[INDEX(i-1, j)] + vel_y[INDEX(i-1, j+1)] + vel_y[INDEX(i, j+1)])*0.25f;
				break;
			case VELY: // vel_y (on the bottom)
				v[0] += (vel_x[INDEX(i, j)] + vel_x[INDEX(i-1, j)] + vel_x[INDEX(i-1, j+1)] + vel_x[INDEX(i, j+1)])*0.25f;
				v[1] += vel_y[INDEX(i, j)];
				break;
			default: // error
				printf("Error in SWSolver::advect\n\twrong FIELD=%d",QUANTITY);
				exit(1);
			}

			// backtrace position
			float srcpi = (float)i - v[0] * dt * 1.0/dx;
			float srcpj = (float)j - v[1] * dt * 1.0/dx;

			// clamp range of accesses (dont access boundaries (at 0 and res-1))
			if(srcpi<0.) srcpi = 1.f;
			if(srcpj<0.) srcpj = 1.f;
			if(srcpi>xRes-1.) srcpi = xRes-2.f;
			if(srcpj>yRes-1.) srcpj = yRes-2.f;

			// interpolate source value
			noi[INDEX(i, j)] = interpolate(array, srcpi, srcpj);
		}
		// grid->switchOldNew(QUANTITY);
}

void SWSolver::updateHeight(){
	int xRes = grid->xRes;
	int yRes = grid->yRes;
	float dx = grid->dx;

	float *eta    = grid->oldFields[ETA];
	float *height = grid->oldFields[HEIGHT];
	float *ground  = grid->oldFields[GROUND];
	const float *vel_x = grid->oldFields[VELX];
	const float *vel_y = grid->oldFields[VELY];

	for (int i = 1; i < xRes - 1; i++)
	for (int j = 1; j < yRes - 1; j++){	
		eta[INDEX(i, j)] += -eta[INDEX(i, j)] * dt * ((vel_x[INDEX(i + 1, j)] - vel_x[INDEX(i, j)]) / dx + (vel_y[INDEX(i, j + 1)] - vel_y[INDEX(i, j)]) / dx);
		if (eta[INDEX(i, j)] < 0.0) 
			eta[INDEX(i, j)] = 0.0;
		height[INDEX(i, j)] = eta[INDEX(i, j)] + ground[INDEX(i, j)];
	}
	
}

void SWSolver::updateVelocity(){
	int xRes = grid->xRes;
	int yRes = grid->yRes;
	float dx = grid->dx;
	float *vel_x = grid->oldFields[VELX];
	float *vel_y = grid->oldFields[VELY];
	float *height = grid->oldFields[HEIGHT];

	for(int i=2; i<xRes-1; i++)
		for(int j=1; j<yRes-1; j++){
			vel_x[INDEX(i, j)] += (g*((height[INDEX(i-1, j)]-height[INDEX(i, j)])/dx)+a_ext[0])*dt;
		}
	for(int i=1; i<xRes-1; i++)
		for(int j=2; j<yRes-1; j++){
			vel_y[INDEX(i, j)] += (g*((height[INDEX(i, j-1)]-height[INDEX(i, j)])/dx)+a_ext[0])*dt;
		}
}

void SWSolver::setBoundary(){
	int xRes = grid->xRes;
	int yRes = grid->yRes;
	float *height = grid->oldFields[HEIGHT];

	for(int i=0; i<xRes; i++){
		height[INDEX(i, 0)]      = height[INDEX(i, 1)];				// lower boundary
		height[INDEX(i, yRes-1)] = height[INDEX(i, yRes-2)];	// upper boundary
	}
	for(int j=0; j<yRes; j++){
		height[INDEX(0, j)]      = height[INDEX(1, j)];				// left boundary
		height[INDEX(xRes-1, j)] = height[INDEX(xRes-2, j)];	// right boundary
	}
}

void SWSolver::advanceTimestep(float _dt){
			dt = _dt;
			advect(ETA);
			grid->switchOldNew(ETA);
			advect(VELX);
			grid->switchOldNew(VELX);
			advect(VELY);
			grid->switchOldNew(VELY);
			updateHeight();
			updateVelocity();
			setBoundary();

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
	float *height = grid->oldFields[HEIGHT];
	float *eta = grid->oldFields[ETA];
	float *ground = grid->oldFields[GROUND];
	float xRes = grid->xRes;
	float yRes = grid->yRes;
	for(int i=0; i<xRes; i++)
		for(int j=0; j<yRes; j++){
			height[INDEX(i, j)] = eta[INDEX(i, j)] + ground[INDEX(i, j)];
		}
}

float SWSolver::computeKineticEnergy(){
	float *vel_x = grid->oldFields[VELX];
	float *vel_y = grid->oldFields[VELY];
	float energy = 0.;
	for (int i = 0; i < grid->xRes; ++i)
	for (int j = 0; j < grid->yRes; ++j)
	{
		float u = vel_x[i+j*grid->xRes];
		float v = vel_y[i+j*grid->xRes];
		energy += 0.5 * (u*u + v*v);
	}
	return energy;
}

float SWSolver::computePotentialEnergy(){
	float *height = grid->oldFields[HEIGHT];
	float energy = 0.;
	for (int i = 0; i < grid->xRes; ++i)
	for (int j = 0; j < grid->yRes; ++j){
		energy += g*height[i+j*grid->xRes];
	}
	return energy;
}


SWRBSolver::SWRBSolver(Sw_grid *_grid, Box *b) : 
	box(b),
	rbs(b),
    grid(_grid),
	alpha(1.0f),
	rho(100.0f),
	g(9.81f)

{
	int xRes = grid->xRes;
	int yRes = grid->yRes;
	displ_old.resize(xRes*yRes);
	displ_new.resize(xRes*yRes);
}

Box* SWRBSolver::getBody(){
	return this->box;
}

void SWRBSolver::advanceTimestep(float dt){
	float *height = grid->oldFields[HEIGHT];
	float *vel_x = grid->oldFields[VELX];
	float *vel_y = grid->oldFields[VELY];
	float dx = grid->dx;
	cout << "setting displacements..." << endl;
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
	bool* isIntersecting = new bool[(x_max-x_min)*(y_max-y_min)];
	std::vector<Vector3f> positions;
	Vector3f r;
	// set the displacements
	int index = 0;
	float totalDisplacement = 0.0f;
	for(int i=x_min; i<x_max; i++)
		for(int j=y_min; j<y_max; j++){
			isIntersecting[index] = calculateDisplacement(i, j, displ_new[INDEX(i, j)], r);
			totalDisplacement += displ_new[INDEX(i, j)];
			if(isIntersecting[index])
				positions.push_back(r);
			index++;
		}

//	cout << "handling body interaction..." << endl;
	cout << totalDisplacement << endl;
	std::vector<Vector3f> forces;
	index = 0;
	for(int i=x_min; i<x_max; i++)
		for(int j=y_min; j<y_max; j++){
			// handle body --> water
			// set height of neighbouring cells
			float dh = (displ_new[INDEX(i, j)]-displ_old[INDEX(i, j)]);
			//height[INDEX(i, j)] -= dh*alpha;
			height[INDEX(i+1, j)] += dh*0.25*alpha;
			height[INDEX(i, j+1)] += dh*0.25*alpha;
			height[INDEX(i-1, j)] += dh*0.25*alpha;
			height[INDEX(i, j-1)] += dh*0.25*alpha;

			// handle water --> body
			if(isIntersecting[index]){
				float hb = 1;  // what is hb in Thue07 eq (11)? 
				float buoyancy = g*dx*dx*displ_new[INDEX(i, j)]*rho;
				forces.push_back(Vector3f(hb*vel_x[INDEX(i, j)], hb*vel_y[INDEX(i, j)], buoyancy));
			}
			index++;
		}

	// add gravity
	forces.push_back(Vector3f(0.0f, 0.0f, -1.0f));
	positions.push_back(box->x);

/*
	// some funny rotation
	forces.push_back(Vector3f(0.0f, 0.0f, 0.0001f));
	positions.push_back(box->x-box->x0/2.0f-box->y0/2.0f-box->z0/2.0f);
*/

	rbs.advanceTimestep(dt, forces, positions);
	forces.clear();
	positions.clear();

	displ_old = displ_new; // note that the vector class has an overloaded "=" which copies the content
}

void SWRBSolver::estimateIndices(Vector3f vertices[8], int &x_min, int &x_max, int &y_min, int &y_max){
	int res[] = {
		grid->xRes,
		grid->yRes
	};
//	float *height = grid->oldFields[HEIGHT];
//	std::vector<Vector3f> convHull;

	// for more accurate estimate of indices
/*	convHull = getConvexHull8XY(vertices);
	for(unsigned int i=0; i<convHull.size(); i++)
		cout << convHull[i][0] << " " << convHull[i][1] << " " << convHull[i][2] << endl;
*/

	// calculate bounding box on grid for a rough estimate of grid points
	// get points with x_max and x_min
	bubbleSortVert(0, vertices, 8);
	x_min = int(vertices[0][0]/grid->dx);		// always round down for minimal index
	x_max = int(vertices[7][0]/grid->dx+1);	// always round up for maximal index
	bubbleSortVert(1, vertices, 8);
	y_min = int(vertices[0][1]/grid->dx);
	y_max = int(vertices[7][1]/grid->dx+1);

	// clamp range of access to boundary values
	if(x_min<0) x_min = 1;
	if(y_min<0) y_min = 1;
	if(x_max>res[0]) x_max = res[0]-2;
	if(y_max>res[1]) y_max = res[1]-2;
}



//! displ is the "height" of the displacement and r is the position on the surface of the body where the intersection happens
// r is later on used to calculate the torque excerted on the body
bool SWRBSolver::calculateDisplacement(int i, int j, float &displ, Vector3f &r){
	int res[] = {
		grid->xRes,
		grid->yRes
	};
	float *height = grid->oldFields[HEIGHT];
	Vector3f P_line(float(i)*grid->dx, float(j)*grid->dx, 0.0f);
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
	for(int k=0; k<3; k++){
		planePoint[k] = vertex_bl;
		planeVector1[k+3] = -planeVector1[k];
		planeVector2[k+3] = -planeVector2[k];
		planePoint[k+3] = vertex_tr;
	}

	int counter = 0;
	int l = 0;
	bool isIntersecting;
	Vector3f P[2];
	// loop through planes and look for intersections
	while(l<6 && counter<2){								// remember that the box is convex, if we change to non-convex boats, counter can be >2
		Vector3f N = (planeVector1[l].cross(planeVector2[l])).normalized();
		if(fabs((V*N)) < 0.001) {							// do not consider the case where the plane is almost perpendicular to the x-y-plane
			l++;											// this will most probably lead to an intersection far away from the rectangle of interest
			continue;
		}
		float t = (planePoint[l]-P_line)*N/(V*N);
		P[counter] = P_line + t*V;

		// if the intersection point is further away from the chosen vertex than p1+p2 its not on the box
		float projectP1 = (P[counter]-planePoint[l])*(planeVector1[l].normalized());
		float projectP2 = (P[counter]-planePoint[l])*(planeVector2[l].normalized());
		float epsilon = -0.001f;
		if(projectP1 < (planeVector1[l].length() + epsilon) &&			
		   projectP2 < (planeVector2[l].length() + epsilon) &&			// is P within a rectangle spanned by |P1| and |P2|?
		   projectP1 > 0.0f &&
		   projectP2 > 0.0f) {											// is P in the "direction" of P1 and P2
			counter++;
		}
		l++;
	}
	// case of intersection (note that the box is convex)
	// for convenience set P[0] to be the lower end (lower z) of the body and P[1] to be the upper end (P[0]<P[1] always)
	if(fabs(P[0][2]-P[1][2]) < 0.0001) counter = 0;	// neglect edges
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
	// case of no intersection or intersecting close to the edges of the box
	else if(counter == 0 || counter == 1){
		displ = 0.0f;
		isIntersecting = false;
	}
	// something weird happened!
	else {
		cout << "weird behaviour" << endl;
        isIntersecting = false;
    }

	return isIntersecting;
}

std::vector<float>* SWRBSolver::getDisplacement(){
	return &this->displ_new;
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
