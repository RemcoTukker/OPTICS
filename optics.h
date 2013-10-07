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

/**
 * The OPTICS algorithm
 */
class Optics {
public:
	typedef float value_t;
	typedef Eigen::Matrix<value_t,1,Eigen::Dynamic> vector_t; // row_vector
	//typedef Eigen::RowVectorXf vector_t;

    //struct to store incomping samples in
	struct DataPoint {
		DataPoint(vector_t sample, size_t ground_truth) {
			this->ground_truth = ground_truth;
			this->data = sample;
		}

        size_t ground_truth;  ///be careful not to touch this until evaluation! (and dont use evaluation results for anything!)
		vector_t data;

        bool processed; // for OPTICS to keep track of which datapoints were processed already

        float reachabilityDistance; // values calculated in OPTICS algorithm
        float coreDistance;         //

	};

	// Cluster keeps track of the samples that belong to one cluster as a result of the algorithm
	struct Cluster {
		std::vector<size_t> r_data;
	};

    //constructor
	Optics(float eps, int minPoints ) {
        minPts = minPoints;
        epsilon = eps;
        dimensionality = -1;
    }

    //destructor doesnt have to do anything
	~Optics() {}

    //this function is executed every once in a while
	void tick() {

        // make sure we make a clean start
        clearResults();

        //execute the OPTICS algorithm
	    startOptics();

	    std::ofstream f;
		std::string file = "./output.data";
		f.open(file.c_str());
		for (int i = 0; i < ordered_set.size(); i++) {
			f << i << " " << ordered_set[i].reachabilityDistance << " " << ordered_set[i].ground_truth << std::endl;
		}
		f.close();

		///TODO: evaluate whether the result is reasonable; otherwise retry with better epsilon and minpts settings!
		//calculate average density and compare it with epsilon (see formula in paper)
		//minPts in reasonable range (10-20)
		//not too much at -1 (undefined / outliers / starts of new clusters) minPts too large
		//not too jagged (invent some measure) minPts too small

	    ///TODO: extract clusters from the results
		//vector setOfSteepDownAreas
		//vector setOfClusters
		float xi = 0.05;
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

		}


	    ///TODO: evaluate result using RAND

	}

	void clearResults() {
        for (int i = 0; i < data_set.size(); i++) {
            data_set[i].processed = false;
            data_set[i].reachabilityDistance = -1;  //any value below zero will do, as distances are always positive
            data_set[i].coreDistance = -1;          // negative values correspond to "undefined"
	    }
        ordered_set.clear();
        clusters.clear();

	}

    //finds the first sample that wasnt processed yet
	void startOptics() {
        for (int i = 0; i < data_set.size(); i++) {
            if (!data_set[i].processed)
                expandClusterOrder(i); //try to expand it into a cluster
	    }
	}

    //takes the number of a sample and attempts to expand it to a cluster
    void expandClusterOrder(const int i) {

        //find neighbouring samples
        std::vector<int> neighbors;
        std::vector<float> distances;
        for (int j = 0; j < data_set.size(); j++) {
            if (j == i) continue; //dont put yourself on your neighbors list..
            vector_t difference = data_set[i].data - data_set[j].data;
            if (difference.norm() < epsilon) {
                neighbors.push_back(j);
                distances.push_back(difference.norm());
            }
        }

        //dont return to this sample anymore
        data_set[i].processed = true;

        //data_set[i].reachabilityDistance stays -1

        //calculate core distance
        if (neighbors.size() >= minPts) {
            std::sort(distances.begin(), distances.end());
            data_set[i].coreDistance = distances[minPts - 1];
        }
        //else it stays -1, undefined

        ordered_set.push_back(data_set[i]); //insert data point in ordered list: start of a new cluster (if its a core object)
        if (data_set[i].coreDistance < 0) return; // if this point is not a core object, return and try again with a new starting point

        //vector with the seeds for expanding the cluster
        std::vector<int> seeds;

        //update the order of the seeds for expanding the cluster
        updateOrderSeeds(neighbors, i, seeds);

        while (seeds.size() > 0) {  ///TODO: try to restructure this, would be much prettier as a recursive function!
            //find the seed with the smallest reachability distance from the seeds vector
            float smallestReachabilityDistance = data_set[seeds[0]].reachabilityDistance;
            int nextSeed = 0;
            for (int j = 1; j < seeds.size(); j++) {
                if (data_set[seeds[j]].reachabilityDistance < smallestReachabilityDistance) {
                    nextSeed = j;
                    smallestReachabilityDistance = data_set[seeds[j]].reachabilityDistance;
                }
            }

            //find neighbours of seed //data_set[seeds[0]];
            std::vector<int> seedsNeighbors;
            std::vector<float> seedsDistances;
            for (int j = 0; j < seeds.size(); j++) {
                if (j == nextSeed) continue;
                vector_t difference = data_set[seeds[j]].data - data_set[seeds[nextSeed]].data;
                if (difference.norm() < epsilon) {
                    seedsNeighbors.push_back(j);
                    seedsDistances.push_back(difference.norm());
                }
            }

            data_set[seeds[nextSeed]].processed = true;

            if (seedsNeighbors.size() >= minPts) {
                std::sort(seedsDistances.begin(), seedsDistances.end());
                data_set[seeds[nextSeed]].coreDistance = seedsDistances[minPts - 1];
            }

            ordered_set.push_back(data_set[seeds[nextSeed]]);

            if (data_set[seeds[nextSeed]].coreDistance < 0) {
				seeds.erase(seeds.begin() + nextSeed);
            	continue; //this point is not a core object, go on with next sample
            }

            updateOrderSeeds(seedsNeighbors, seeds[nextSeed], seeds); //expand the cluster
            seeds.erase(seeds.begin() + nextSeed);
        }

    }

    void updateOrderSeeds(const std::vector<int> &neighbors, const int i, std::vector<int> &seeds) {

        for (int j = 0; j < neighbors.size(); j++) {
            if (data_set[neighbors[j]].processed) continue;
                //hrmm.. this is only here 'cause we never take stuff out of the seeds list.... ugly!

            vector_t difference = data_set[i].data - data_set[neighbors[j]].data;
            float newReachibilityDistance = std::max(data_set[i].coreDistance, difference.norm());

            //insert neighbour as seed if it isnt in the list yet
            if (data_set[neighbors[j]].reachabilityDistance < 0) {
                seeds.push_back(neighbors[j]);
                data_set[neighbors[j]].reachabilityDistance = newReachibilityDistance;
            } else { //should be in the list already, just check if we have to update r-distance
                data_set[neighbors[j]].reachabilityDistance = std::min(newReachibilityDistance, data_set[neighbors[j]].reachabilityDistance);
            }
        }
    }

    //this function is executed when a new sample is coming in: vector with values, then a label
	void addSample(std::vector<value_t> & x, size_t label) {

        if (dimensionality == -1) {
            dimensionality = x.size();
            std::cout << "First sample recieved, dimensionality of data set to " << dimensionality << std::endl;
        }

        if (x.size() != dimensionality) { ///TODO: handle this more gracefully and just return a warning
            std::cout << "Error: dimensionality of incoming sample incorrect." << std::endl;
            std::cout << "This sample has D = " << x.size() << " while the first sample had D = " << dimensionality << std::endl;
            return;
        }

        //copy the data to an Eigen row vector
		vector_t v(x.size());  //v.setZero();
		v = vector_t::Map(x.data(), 1, x.size()); // create row-vector

        //store it in data_set
		data_set.push_back(DataPoint(v, label));

        ///Preferably process each datapoint as it comes in instead of at the "end" in the tick function!

	}

	///TODO: add overloaded addSample function to accept samples without label
	///        (reserve label 0 for unknown or use some different solution?)


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


protected:

private:

    float epsilon;
    int minPts;
    int dimensionality;

	std::vector<DataPoint> data_set;    //every sample thats coming in is stored here
	std::vector<DataPoint> ordered_set;  //the result of the basic OPTICS algorithm
	std::vector<Cluster> clusters;      //using some heuristics, clusters are made out of that

};


#endif /* OPTICS_H_ */
