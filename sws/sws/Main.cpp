#include "util.h"

const int N_STEPS = 100;

int main(){

	for(int t=0; t<N_STEPS; t++){
		advectHeight(eta, vel);
		advectVelX(vel_x, vel);
		advectVelY(vel_y, vel);

		updateHeight(eta, vel);
		updateVel(height, vel_x, vel_y);
	}
}