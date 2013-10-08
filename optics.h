/**
 *
 * @brief ...
 * @file optics.h
 *
 * This file is created at Almende B.V. and Distributed Organisms B.V. It is open-source software and belongs to a
 * larger suite of software that is meant for research on self-organization principles and multi-agent systems where
 * learning algorithms are an important aspect.
 *
 * This software is published under the GNU Lesser General Public license (LGPL).
 *
 * It is not possible to add usage restrictions to an open-source license. Nevertheless, we personally strongly object
 * against this software being used for military purposes, factory farming, animal experimentation, and "Universal
 * Declaration of Human Rights" violations.
 *
 * Copyright (c) 2013 Remco Tukker <remco@almende.org>
 *
 * @author    R Tukker
 * @date      Oct 6, 2013
 * @project   Replicator
 * @company   Almende B.V.
 * @company   Distributed Organisms B.V.
 * @case      Clustering
 */

#ifndef OPTICS_H_
#define OPTICS_H_

#include "Eigen/LU" // <Eigen/LU>
#include <vector>
#include <algorithm>
#include <fstream>
#include <map>

/**
 * A class implementing the OPTICS algorithm
 *
 * Known issues:
 *  - Can only handle 32k data points for now (as I use ints to keep track of individual data points)
 *  - There's no way to remove (old) data points or send a data point only for asking in which cluster it belongs
 *  - There's currently no way to use a couple of labels from input data for better naming of clusters
 *  - OPTICS is a density-based clustering algorithm; dont expect good results when clusters/classes overlap
 *
 * TODO:
 *  - Finish algorithms
 *  - Fix known issues (as far as possible)
 *  - Fix other TODO's scattered in the code
 *  - Use an indexer to speed up neighborhood queries...
 */

class Optics {
public:
	typedef float value_t;
	typedef Eigen::Matrix<value_t,1,Eigen::Dynamic> vector_t; //typedef Eigen::RowVectorXf vector_t;

    //struct to store incomping samples in
	struct DataPoint {
		DataPoint(vector_t sample, size_t ground_truth) {
			this->ground_truth = ground_truth;
			this->data = sample;
		}

        size_t ground_truth;  ///do not touch this until evaluation! (and dont use evaluation results for anything!)
		vector_t data;

        bool processed;             // for OPTICS to keep track of which datapoints were processed already
        float reachabilityDistance; // values calculated in OPTICS algorithm
        float coreDistance;         //

	};

	///TODO: adapt this to our need...
	// Cluster keeps track of the samples that belong to one cluster as a result of the algorithm
	struct Cluster {
		std::vector<size_t> r_data;
	};

    //constructor
	Optics(float eps, int minPoints, bool fixDimensionMismatches = false ) {
        minPts = minPoints;
        epsilon = eps;
        dimensionality = -1;
        permissive = fixDimensionMismatches;
    }

    //destructor doesnt have to do anything
	~Optics() {}

    //this function is executed every once in a while
	void tick() {

        // make sure we make a clean start
        clearResults();

        //execute the OPTICS algorithm
	    startOptics();

		///TODO: remove this
		//give output of OPTICS as a data file for now
	    std::ofstream f;
		std::string file = "./output.data";
		f.open(file.c_str());
		for (int i = 0; i < ordered_set.size(); i++) {
			f << i << " " << ordered_set[i].reachabilityDistance << " " << ordered_set[i].ground_truth << std::endl;
		}
		f.close();

		//evaluate results to see whether parameters were reasonable
		evaluateResults();

		//extract clusters from the results
		extractClusters();

	    //compare result with ground truth using RAND
	    compareWithGroundTruth();

	}

    //this function is executed when a new sample is coming in: vector with values, then a label
	void addSample(std::vector<value_t> & x, size_t label) {

		///TODO: make chit-chat optional / log / ... (?)

        if (dimensionality == -1) {
            dimensionality = x.size();
            std::cout << "First sample recieved, dimensionality of data set to " << dimensionality << std::endl;
        }

		if (x.size() != dimensionality) {
			if (!permissive) {
				std::cout << "Warning: incoming sample has incorrect dimension! Sample is rejected!" << std::endl;
				std::cout << "This sample has D = " << x.size()
								<< " while the first sample had D = " << dimensionality << std::endl;
				return;
			}

			if (x.size() < dimensionality) {
				std::cout << "Warning: dimensionality of incoming sample too low!" << std::endl;
				std::cout << "This sample has D = " << x.size()
								<< " while the first sample had D = " << dimensionality << std::endl;
				std::cout << "We will assume that the remaining dimensions default to zero..." << std::endl;
				x.resize(dimensionality);
			} else if (x.size() > dimensionality) {
				std::cout << "Warning: dimensionality of incoming sample too high!" << std::endl;
				std::cout << "This sample has D = " << x.size()
								<< " while the first sample had D = " << dimensionality << std::endl;
				std::cout << "We will ignore the extra dimensions..." << std::endl;
				x.resize(dimensionality);
			}
		}


        //copy the data to an Eigen row vector
		vector_t v(x.size());  //v.setZero();
		v = vector_t::Map(x.data(), 1, x.size()); // create row-vector

        //store it in data_set
		data_set.push_back(DataPoint(v, label));

        ///Preferably process each datapoint as it comes in instead of at the "end" in the tick function!
        /// is it possible with OPTICS without running the whole algorithm from start to end?

	}

	///TODO: add overloaded addSample function to accept samples without label
	///        (reserve label 0 or -1 for unknown or use some different solution?)

protected:

private:

    float epsilon;
    int minPts;
    int dimensionality;

    bool permissive;

    std::vector<int> seeds;
	std::vector<DataPoint> data_set;    //every sample thats coming in is stored here
	std::vector<DataPoint> ordered_set;  //the result of the basic OPTICS algorithm
	std::vector<Cluster> clusters;      //using some heuristics, clusters are made out of that

	/**
	 *   Preparing for a new run of the algorithm
	 **/
	void clearResults() {
        for (int i = 0; i < data_set.size(); i++) {
            data_set[i].processed = false;
            data_set[i].reachabilityDistance = -1;  //any value below zero will do, as distances are always positive
            data_set[i].coreDistance = -1;          // negative values correspond to "undefined"
	    }
        ordered_set.clear();
        clusters.clear();
        seeds.clear();

	}

	void evaluateResults() {

		///TODO: evaluate whether the result is reasonable; otherwise retry with better epsilon and minpts settings!
		//calculate average density and compare it with epsilon (see formula in paper) (perhaps do this beforehand)
		//minPts in reasonable range (10-20)
		//not too much at -1 (undefined / outliers / starts of new clusters) minPts too large
		//not too jagged (invent some measure) minPts too small
	}

	void compareWithGroundTruth() {
		///TODO: use RAND for comparison
		///possibly dont have this function here at all, its not the right place (?)
		///only for testing / comparison purposes I guess..
			/*
	//get an evaluation of the result, given a ground truth..
	void evaluate() {
		for (int k = 0; k < clusters.size(); ++k) {
			for (int i = 0; i < clusters[k].r_data.size(); ++i) {
				labels[clusters[k].r_data[i]].prediction = k;
//				std::cout << "Cluster " << k << ", data item " << clusters[k].r_data[i] << std::endl;
			}
		}

		size_t a, b, c, d; a = b = c = d = 0;
		for (int i = 0; i < data_set.size(); ++i) {
			for (int j = i+1; j < data_set.size(); ++j) {
				if (labels[i].ground_truth == labels[j].ground_truth) {
					if (labels[i].prediction == labels[j].prediction) {
						a++;
					} else {
						c++;
					}
				} else {
					if (labels[i].prediction == labels[j].prediction) {
						d++;
					} else {
						b++;
					}
				}
			}
		}

		value_t quality = (a+b) / value_t(a+b+c+d);
		//std::cout << "Rand index is (" << a << "+" << b << ") /" << a+b+c+d << " = " << quality << std::endl;
		std::cout << "Rand index is " << quality << std::endl;
	}
	*/


	}

	void extractClusters() {

		///TODO: extract clusters from the results -> use clusters vector
		///TODO: perhaps add a DBSCAN mode for extracting DBSCAN results
		//vector setOfSteepDownAreas
		//vector setOfClusters
		/*float xi = 0.05;
		int index = 0;
		float mib = 0;
		while (index < ordered_set.size() - 1) {
			mib = std::max(mib, ordered_set[index].reachabilityDistance);

			bool startDown = false;
			bool startUp = false;

			//check if we have to do with a steep point:
			if (ordered_set[index].reachabilityDistance < 0 && ordered_set[index + 1].reachabilityDistance < 0) {
				//two undefined points; current one is an outlier
			} else if (ordered_set[index].reachabilityDistance < 0 && ordered_set[index+1].reachabilityDistance >= 0) {
				//going from undefined to a distance; definitely start of a steep down area!
				startDown = true;
			} else if (ordered_set[index].reachabilityDistance >= 0 && ordered_set[index+1].reachabilityDistance < 0) {
				//going from distance to undefined point; definitely start (and end) of steep up area
				startUp = true;
			} else if (ordered_set[index].reachabilityDistance *(1-xi) >= ordered_set[index+1].reachabilityDistance) {
				//start of steep down area!
				startDown = true;
			} else if (ordered_set[index].reachabilityDistance <= ordered_set[index+1].reachabilityDistance *(1-xi)) {
				//start of steep up area
				startUp = true;
			} else {
				//nothing special.. point of a cluster

			}


			if (startDown) //start of steep down area D at index
			{
				//update mib-values and filter SetOfSteepDownAreas(*)
				//set D.mib = 0
				//add D to the SetOfSteepDownAreas


				//index = end of D + 1; mib = r(index)


			}
			else if (startUp) //start of steep up area U at index
			{
				//update mib-values and filter SetOfSteepDownAreas


				//index = end of U + 1; mib = r(index)


				//FOR EACH D in SetOfSteepDownAreas DO
				//	IF(combination of D and U is valid AND(**)
				//		satisfies cluster conditions 1, 2, 3a)
				//		compute [s, e] add cluster to SetOfClusters

			}
			else {
				index++;
			}

		}*/
	}



	/****************************************************
	 * Function to start OPTICS and keep expanding seeds
	 ****************************************************/
	void startOptics() {
        for (int i = 0; i < data_set.size(); i++) {
            if (!data_set[i].processed) {
            	seeds.push_back(i); //finds the first sample that wasnt processed yet
				expandSeeds();      //and then expands it to a cluster, if possible
			}
	    }
	}

	/**
	 * Helper function to find the seed that should be processed next: the one with the smallest reachability distance
	 */

    int seedWithSmallestRDistance() {

    	float smallestReachabilityDistance = data_set[seeds[0]].reachabilityDistance;
    	int seed = 0;

		for (int j = 1; j < seeds.size(); j++) {
            if ( data_set[seeds[j]].reachabilityDistance < 0 ) continue; //undefined, so not smaller than anything
            if ( (smallestReachabilityDistance < 0) //anything is smaller than undefined
					|| (data_set[seeds[j]].reachabilityDistance < smallestReachabilityDistance) ) {
				seed = j;
				smallestReachabilityDistance = data_set[seeds[j]].reachabilityDistance;
            }
        }

		return seed;
    }

    /***********************************************************
     * This function is the actual OPTICS algorithm
     * Make sure you have at least 1 seed before calling it
     ***********************************************************/
    void expandSeeds() {

        //first find seed thats up next; thats the one with the smallest reachability distance
        int nextSeed = seedWithSmallestRDistance();
        int i = seeds[nextSeed];  //store the index of the sample we are interested in
        seeds.erase(seeds.begin() + nextSeed); //and take it off the seeds list

        //then find the seeds neighbors
        std::map<int, float> neighbors;
        std::vector<float> distances;
        for (int j = 0; j < data_set.size(); j++) {  /// TODO: this is the main/only candidate for speedup by indexing
            if (j == i) continue;                    //  dont put yourself on your neighbors list..
            float distance = (data_set[i].data - data_set[j].data).norm();
            if (distance >= epsilon) continue;       //  point is too far away to be a neighbor
            neighbors.insert( std::pair<int, float>(j, distance) );
            distances.push_back(distance);           // this one only for sorting in the next step
        }

        //calculate core distance of current seed (it stays -1 if its not a core object)
        if (neighbors.size() >= minPts) {
            std::partial_sort(distances.begin(), distances.begin() + minPts ,distances.end());
            data_set[i].coreDistance = distances[minPts - 1];
        }

        //do some bookkeeping...
        data_set[i].processed = true;             //dont return to this sample anymore
        ///TODO: better to only store order with ints instead of whole structs
        ordered_set.push_back(data_set[i]);       //insert data point in ordered list
        if (data_set[i].coreDistance < 0) return; //if this point is not a core object, return and try next seed

        //finally, we have a core object, so add its neighbours to the seed list
        for (std::map<int, float>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
            if (data_set[it->first].processed) continue; //dont add neighbors that have been processed already

            //check if we should update reachability distance from current sample
            float newReachibilityDistance = std::max(data_set[i].coreDistance, it->second);

            //insert neighbour as seed if it isnt in the list yet
            if (data_set[it->first].reachabilityDistance < 0) {
                seeds.push_back(it->first);
                data_set[it->first].reachabilityDistance = newReachibilityDistance;
            } else { //should be in the list already, just check if we have to update r-distance
                data_set[it->first].reachabilityDistance =
						std::min(newReachibilityDistance, data_set[it->first].reachabilityDistance);
            }
        }

        //keep going until there are no seeds left
        if (!seeds.empty()) expandSeeds();
    }

};


#endif /* OPTICS_H_ */
