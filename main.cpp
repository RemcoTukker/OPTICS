#include <iostream>
#include "optics.h"
#include <iostream>
#include <fstream>
#include "data.hpp"

//using namespace std;

int main()
{

	//load data

	typedef float value_t;
	data<value_t> d;
	int index;



	std::ifstream f;
	//std::string file = "./data/gauss.data";
	std::string file = "./data/abalone3.data";
	//std::string file = "./data/iris2.data";

// obtained # of clusters by $(cat ../../data/abalone3.data | cut -f11 -d',' | sort -n | uniq | wc -l)
//predefined_clusters = 28;
//file = "../../data/gaussian3d1.data";
//predefined_clusters = 3;
//		file = "../../data/gaussian1.data";
//		predefined_clusters = 2;

	f.open(file.c_str());
	if (f.is_open()) {
		f >> d;
	} else {
		std::cerr << "File " << file << " does not exist " << std::endl;
	}

	index = 0;

	//feed stuff into optics
	std::vector<value_t> & item = d.pop();
	int D = item.size() - 1; // dimensionality of data samples (length of vector minus field for label)
	std::cout << "Dimensionality is " << D << std::endl;

	//int K = predefined_clusters;
	float eps = 0.7;
	int minPts = 25;
	Optics optics(eps, minPts);

	int S = d.size(); // # samples
	std::cout << "Load all " << S << " samples" << std::endl;

	for (int s = 0; s < S; ++s) {

		int label = item[D];
		item.pop_back(); //remove label from data vector
		optics.addSample(item, label);
		item = d.pop();
	}

	optics.tick();

	std::cout << std::endl;

	//kmeans.evaluate();




    std::cout << "Hello world!" << std::endl;
    return 0;
}
