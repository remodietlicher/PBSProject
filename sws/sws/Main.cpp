#include "util.h"
#include "sws.h"
#include "viewer.h"

#include <fstream>
#include <iostream>
#include <math.h>

const int N_STEPS = 1;

using namespace std;

int main(int argc, char** argv){
	int res[2] = {100, 100};
	SWSolver sws(res[0], res[1], 1, 1, 0.001f);
	SWViewer swv(&sws);
	std::vector<std::vector<float>> heightmaps;
	std::vector<float> initEta;
	initEta.resize(res[0]*res[1]);
	for(int i=0; i<res[0]; i++)
		for(int j=0; j<res[1]; j++){
		if(i<4*res[0]/6 && i>2*res[0]/6 && j<4*res[1]/6 && j>2*res[1]/6) initEta[INDEX(i, j)] = 1.1f;
		else initEta[INDEX(i, j)] = 1.0f;
	}

	sws.setEta(initEta);

	sws.advanceTimestep();
	swv.drawScene();


/*	ofstream data_heightmap;
	data_heightmap.open("testdata//heightmap.txt");
	int l=0;
	int a = 1;
	while(a == 1){
		sws.advanceTimestep();

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
	*/
}