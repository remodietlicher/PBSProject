#ifndef UTIL_H
#define UTIL_H

#include "assert.h"
#include <math.h>

class Vec3f{
private:
	float xyz[3];
public:
	Vec3f() {
		this->xyz[0] = 0.0f;
		this->xyz[1] = 0.0f;
		this->xyz[2] = 0.0f;
	}

	Vec3f(const float x, const float y, const float z){
		this->xyz[0] = x;
		this->xyz[1] = y;
		this->xyz[2] = z;
	}
	float& operator() (int i){ 	return xyz[i]; }
	Vec3f operator*(float s){ 	return Vec3f(xyz[0]*s, xyz[1]*s, xyz[2]*s); }
	Vec3f operator/ (float s){ 	return Vec3f(xyz[0]/s, xyz[1]/s, xyz[2]/s); }
	Vec3f operator+ (Vec3f v){ 	return Vec3f(xyz[0] + v(0), xyz[1] + v(1), xyz[2] + v(2)); }
	Vec3f operator- (Vec3f v){ 	return Vec3f(xyz[0] - v(0), xyz[1] - v(1), xyz[2] - v(2)); }
	float operator* (Vec3f v){ 	return xyz[0]*v(0) + xyz[1]*v(1) + xyz[2]*v(2); }
	Vec3f cross(Vec3f v){ 	return Vec3f(xyz[2]*v(3) - xyz[3]*v(2), xyz[3]*v(1) - xyz[1]*v(3), xyz[1]*v(2) - xyz[2]*v(1)); }
	float length(){ return sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]); }
	Vec3f normalize(){ return *this/this->length(); }
};

#endif