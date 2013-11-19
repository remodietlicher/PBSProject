#include "util.h"
#include "sws.h"

#include <fstream>
#include <iostream>
#include <math.h>

const int N_STEPS = 1000;

using namespace std;

int main(){
	int xRes = 10;
	int yRes = 10;
	SWSolver sws(xRes, yRes, 1, 1, 0.01f);
	std::vector<std::vector<float>> heightmap;
	std::vector<float> initEta;
	initEta.resize(xRes*yRes);
	for(int i=0; i<xRes*yRes; i++){
		initEta[i] = (float)i/(float)(xRes*yRes);
	}

	sws.setEta(initEta);

	for(int t=0; t<N_STEPS; t++){
		sws.advect(ETA);
		sws.advect(VELOCITY_X);
		sws.advect(VELOCITY_Y);
		sws.updateHeight();
		sws.updateVelocity();
		heightmap.push_back(sws.getHeightMap());
	}

	ofstream data_heightmap;
	data_heightmap.open("testdata//heightmap.txt");
	for(int i=0; i<xRes; i++)
		for(int j=0; j<yRes; j++){
			data_heightmap << heightmap[N_STEPS-1][INDEX(i, j)] << endl;
		}
	data_heightmap.flush();
}