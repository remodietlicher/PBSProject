#include "util.h"
#include "sws.h"

#include <fstream>
#include <iostream>
#include <math.h>

const int N_STEPS = 1;

using namespace std;

int main(){
	int res[2] = {20, 20};
	SWSolver sws(res[0], res[1], 1, 1, 0.01f);
	std::vector<std::vector<float>> heightmaps;
	std::vector<float> initEta;
	initEta.resize(res[0]*res[1]);
	for(int i=0; i<res[0]; i++)
		for(int j=0; j<res[1]; j++){
		if(i<4*res[0]/6 && i>2*res[0]/6 && j<4*res[1]/6 && j>2*res[1]/6) initEta[INDEX(i, j)] = 1.1f;
		else initEta[INDEX(i, j)] = 1.0f;
	}

	sws.setEta(initEta);

	ofstream data_heightmap;
	data_heightmap.open("testdata//heightmap.txt");
	int l=0;
	int a = 1;
	while(a == 1){
		cout << "advecting eta..." << endl;
		sws.advect(ETA);
		cout << "advecting velocity_x..." << endl;
		sws.advect(VELOCITY_X);
		cout << "advecting velocity_y..." << endl;
		sws.advect(VELOCITY_Y);

		cout << "updating heights..." << endl;
		sws.updateHeight();

		cout << "updating velocities..." << endl;
		sws.updateVelocity();

		cout << "setting boundaries..." << endl;
		sws.setBoundary();

		cout << "writing heightmap to file..." << endl;
		heightmaps.push_back(sws.getHeightMap());
		for(int i=0; i<res[0]; i++)
			for(int j=0; j<res[1]; j++){
				data_heightmap << heightmaps[l][INDEX(i, j)] << endl;
			}
		data_heightmap.flush();
		l++;
		cout << "continue? (1/0) ";
		cin >> a;
	}
}