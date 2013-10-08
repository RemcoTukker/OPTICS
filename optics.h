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

	// Cluster stores a cluster result as indices to ordered_set
	struct Cluster {
		int start;
		int end;
	};

    //constructor
	Optics(float eps, int minPoints, bool fixDimensionMismatches = false ) {

        ///TODO: check
        minPts = minPoints;
        epsilon = eps;

        dimensionality = -1;
        permissive = fixDimensionMismatches;
    }

    //destructor doesnt have to do anything
	~Optics() {}

    //this function is executed every once in a while
	void tick() {

        clearResults(); // make sure we make a clean start
	    startOptics(); //execute the OPTICS algorithm

		///TODO: remove this when not necessary anymore
		//give output of OPTICS as a data file for now
	    std::ofstream f;
		std::string file = "./output.data";
		f.open(file.c_str());
		for (int i = 0; i < ordered_set.size(); i++) {
			f << i << " " << ordered_set[i].reachabilityDistance << " " << ordered_set[i].ground_truth << std::endl;
		}
		f.close();

		evaluateResults(); //evaluate results to see whether parameters were reasonable
		extractClusters(); //extract clusters from the results

		for (int i = 0; i < clusters.size(); i++) {
			std::cout << i << " " << clusters[i].start << " " << clusters[i].end << std::endl;
		}

	    compareWithGroundTruth(); //compare result with ground truth using RAND

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

	enum SteepAreas { STEEPUP, STEEPDOWN, CLUSTER, OUTLIER };

	SteepAreas checkPoint(int i, float xi) {
		if ((ordered_set[i].reachabilityDistance < 0) && (ordered_set[i+1].reachabilityDistance < 0)) return OUTLIER;
		if ((ordered_set[i].reachabilityDistance >= 0) && (ordered_set[i+1].reachabilityDistance < 0)) return STEEPUP;
		if (ordered_set[i].reachabilityDistance <= ordered_set[i+1].reachabilityDistance *(1-xi)) return STEEPUP;
		if ((ordered_set[i].reachabilityDistance < 0) && (ordered_set[i+1].reachabilityDistance >= 0)) return STEEPDOWN;
		if (ordered_set[i].reachabilityDistance *(1-xi) >= ordered_set[i+1].reachabilityDistance) return STEEPDOWN;
		return CLUSTER;
	}

	int findEndOfSteepRegion(int i, SteepAreas updown, float xi) {
		int endPoint = i;
		int flatPointsCounter = 0;
		for (int j = i + 1; j < ordered_set.size() - 1; j++) {

			//going in reverse direction stops steep area
			if ((updown==STEEPUP)&&(ordered_set[j].reachabilityDistance < ordered_set[j-1].reachabilityDistance)) break;
			if ((updown==STEEPDOWN)&&(ordered_set[j].reachabilityDistance>ordered_set[j-1].reachabilityDistance)) break;

			//check if we have to do with a steep point or a flat point
			if (checkPoint(j, xi) == updown) {
				endPoint = j;
				flatPointsCounter = 0;
			} else {
				flatPointsCounter++;
				if (flatPointsCounter > minPts) break;
			}
		}
		return endPoint;
	}

	struct steepDownArea {
		int start;
		int end;
		float Dmib;
	};

	void tryToMakeACluster(int startUpArea, int endUpArea, std::vector<steepDownArea> &setOfSteepDownAreas, float xi) {

		for (int i = 0; i < setOfSteepDownAreas.size(); i++) {
			//check if we have a valid cluster (condition sc2* in the paper)
			if (setOfSteepDownAreas[i].Dmib > ordered_set[endUpArea].reachabilityDistance * (1 - xi)) continue;

			//determine start and end point
			float reachStart = ordered_set[setOfSteepDownAreas[i].start].reachabilityDistance;
			float reachEnd = ordered_set[endUpArea + 1].reachabilityDistance;
			int SoC = 0;
			int EoC = 0;
			if ( (reachStart * (1 - xi) < reachEnd) && (reachEnd * (1 - xi) < reachStart) ) { //max xi apart
				SoC = reachStart;
				EoC = reachEnd;
			} else if ((reachStart < 0) || (reachStart > reachEnd) ) { //more than xi apart
				EoC = reachEnd;
				float minDistance = fabs(reachEnd - reachStart);
				int closestPoint = setOfSteepDownAreas[i].start;
				for (int j = setOfSteepDownAreas[i].start+1; j <= setOfSteepDownAreas[i].end; j++) { //<= correct here?
					if ( fabs(reachEnd - ordered_set[j].reachabilityDistance) < minDistance ) {
						minDistance = fabs(reachEnd - ordered_set[j].reachabilityDistance);
						closestPoint = j;
					}
				}
				SoC = closestPoint;
			} else {  // (reachStart < reachEnd ) and more than xi apart
				SoC = reachStart;
				float minDistance = fabs(reachEnd - reachStart);
				int closestPoint = startUpArea;
				for (int j = startUpArea+1; j <= endUpArea; j++) { //<= correct here?
					if ( fabs(reachStart - ordered_set[j].reachabilityDistance) < minDistance ) {
						minDistance = fabs(reachStart - ordered_set[j].reachabilityDistance);
						closestPoint = j;
					}
				}
				EoC = closestPoint;
			}

			if (EoC - SoC < minPts) continue; //cluster isnt large enough

			// add cluster to SetOfClusters
			clusters.push_back({SoC, EoC});
		}



	}

	void extractClusters() {

		///TODO: extract clusters from the results -> use clusters vector
		///TODO: perhaps add a DBSCAN mode for extracting DBSCAN results
		std::vector<steepDownArea> setOfSteepDownAreas;
		//vector setOfClusters

		//add a temporary sample at the end with reachability -1
		DataPoint terminator(Eigen::RowVector2f(0, 0), 0);
		terminator.reachabilityDistance = -1;
		ordered_set.push_back(terminator);

		float xi = 0.05;
		int index = 0;
		float mib = 0;
		while (index < ordered_set.size() - 1) { //dont consider the terminator sample
			mib = std::max(mib, ordered_set[index].reachabilityDistance);
			int end;

			switch (checkPoint(index, xi)) {
				case STEEPDOWN:
					///TODO: update mib-values and filter SetOfSteepDownAreas(*)
					end = findEndOfSteepRegion(index, STEEPDOWN, xi);
					setOfSteepDownAreas.push_back( {index, end, 0} );
					index = end + 1;
					mib = ordered_set[index].reachabilityDistance;
					break;
				case STEEPUP:
					///TODO: update mib-values and filter SetOfSteepDownAreas
					end = findEndOfSteepRegion(index, STEEPUP, xi);
					tryToMakeACluster(index, end, setOfSteepDownAreas, xi);
					index = end + 1;
					mib = ordered_set[index].reachabilityDistance;
					break;
				default:
					index++;
			}
		}

		ordered_set.pop_back(); //remove terminator sample again
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

	/******************************************************************************************************************
	 * Helper function to find the seed that should be processed next: the one with the smallest reachability distance
	 *****************************************************************************************************************/

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
