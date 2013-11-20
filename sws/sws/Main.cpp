#include "util.h"
#include "sws.h"

#include <fstream>
#include <iostream>
#include <math.h>

const int N_STEPS = 10;

using namespace std;

int main(){
	int res[2] = {10, 10};
	SWSolver sws(res[0], res[1], 1, 1, 0.05f);
	std::vector<std::vector<float>> heightmaps;
	std::vector<float> initEta;
	initEta.resize(res[0]*res[1]);
	for(int i=0; i<res[0]; i++)
		for(int j=0; j<res[1]; j++){
		if(i<4*res[0]/6 && i>2*res[0]/6 && j<4*res[1]/6 && j>2*res[1]/6) initEta[INDEX(i, j)] = 2.0f;
		else initEta[INDEX(i, j)] = 1.0f;
	}

	sws.setEta(initEta);

	for(int t=0; t<N_STEPS; t++){
		sws.advect(ETA);
		sws.advect(VELOCITY_X);
		sws.advect(VELOCITY_Y);

		float sumBefore = 0;
		for(int i=0; i<res[0]*res[1]; i++)
			sumBefore += sws.getHeightMap()[i];

		sws.updateHeight();

		float sumAfter = 0;
		for(int i=0; i<res[0]*res[1]; i++)
			sumAfter += sws.getHeightMap()[i];

		cout << sumBefore << " " << sumAfter << endl;

		sws.updateVelocity();

		sws.setBoundary();
		heightmaps.push_back(sws.getHeightMap());
	}

	ofstream data_heightmap;
	data_heightmap.open("testdata//heightmap.txt");
	for(int i=0; i<res[0]; i++)
		for(int j=0; j<res[1]; j++){
			data_heightmap << heightmaps[N_STEPS-1][INDEX(i, j)] << endl;
		}
	data_heightmap.flush();
}