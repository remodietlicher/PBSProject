#include "viewer.h"
#include <stdlib.h>

SWViewer::SWViewer(SWSolver *sws){
	this->sws = sws;
	this->normals = new Vec3f[sws->getXRes()*sws->getYRes()];
}

void SWViewer::cleanup(){
	delete normals;
}

void SWViewer::computeNormals(){
	const std::vector<float>& h = sws->getHeightMap();

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
}
int SWViewer::width() 
{ 
	return sws->getXRes(); 
}

int SWViewer::length() 
{ 
	return sws->getYRes(); 
}

Vec3f SWViewer::getNormal(int i, int j)
{
	return normals[V_INDEX(i, j)];
}

float SWViewer::getHeight(int i, int j)
{
	const std::vector<float>& h = sws->getHeightMap();
	return h[V_INDEX(j, i)];
}



