#pragma once
#include <iostream>
#include <string>
#include <vector>
#include "City.h"
#include "Trip.h"
#include <algorithm>
#include <chrono>
#include <numeric>
#include <time.h>
#include <fstream>
#include <unordered_set>
#include <thread>
#include <functional> // for std::ref
#include <map>
#include <set>
#include <math.h>

using namespace std;

bool comparePaths(Trip i1, Trip i2) {
	return(i1.getPathLength() < i2.getPathLength());
}

bool compareMeanDevs(Trip i1, Trip i2) {
	return(i1.getLRSProb() > i2.getLRSProb());
}

float genRandom() { //generates random number between 0 and 1
	return ((float)rand()) / RAND_MAX;
}

double genRandomDouble() { //generates random number between 0 and 1
	return ((double)rand()) / RAND_MAX;
}

//swapping two random cities
void mutate(Trip& gene,int numCities) { 
	int firstSwapIndex = rand() % (numCities- 1);
	int secondSwapIndex = rand() % (numCities - 1);
	vector<City> mutatedPath = gene.getPath();
	City temp = mutatedPath[firstSwapIndex];
	mutatedPath[firstSwapIndex] = mutatedPath[secondSwapIndex];
	mutatedPath[secondSwapIndex] = temp;
	gene.setPath(mutatedPath);
	gene.calcPathLength();
}

//Takes a length, randomizes a selection in the path of that length, scrambles them, then puts them back in the path
void scrambleMutate(Trip& gene,int mutateLength, int numCities) {
	vector<City> path = gene.getPath();
	int startIndex = rand() % (numCities - mutateLength - 1);
	vector<City> cities;
	for (int i = startIndex; i < (startIndex + mutateLength); i++) {
		cities.push_back(path[i]);
	}
	random_shuffle(cities.begin(), cities.end());
	for (int i = startIndex; i < (startIndex + mutateLength); i++) {
		path[i] = cities[i-startIndex];
	}
	gene.setPath(path);
	gene.calcPathLength();

}

/*
Moro mutate works by doing several swaps of cities in a path
These swaps occur at random points, but it involves only swapping two cities 
that are next to each other.  So if we choose index 8, cities 8 and 9 are swapped.
If its 0, cities 0 and 1 are swapped. We do this numSwaps times.
*/
void moroMutate(Trip& gene, int numSwaps,int numCities) {
	vector<City> path = gene.getPath();
	for (int i = 0; i < numSwaps; i++) {
		int startSwap = rand() % (numCities - 1);
		City temp = path[startSwap];
		path[startSwap] = path[startSwap + 1];
		path[startSwap + 1] = temp;
	}
	gene.setPath(path);
	gene.calcPathLength();
}

//Stochastic Universal Sampling
void SUSSelection(vector<Trip>& genePool, vector<Trip>& parents,int popSize) {
	int n = genePool.size();
	double factor = static_cast<double>(1) / n;
	float fHat = 0;
	float a = genRandom();
	
	for (auto& gene : genePool) {
		fHat += gene.getPathLength();
	}

	float mean = factor * fHat;
	float delta = mean * a;
	int j = 0;
	float sum = genePool[j].getPathLength();

	do {
		if (delta < sum) {
			parents.push_back(genePool[j]);
			delta = delta + sum;
		}
		else {
			sum = sum + genePool[j].getPathLength();
			j++;
		}
	} while (j < n);
}

void exponentialRankSelection(vector<Trip>& genePool, vector<Trip>& parents, int selectionPressure) {
	double n = static_cast<double>(genePool.size());
	double c = (n * 2 * (n - 1)) / (6 * (n - 1) + n);
	for (int i = 0; i < n; i++) {
		double p = 1.0f * exp((-i) / c);
		genePool[i].setERSProb(p);
	}

	//not sure about this one
}



//Rank based selection
//https://en.wikipedia.org/wiki/Selection_(genetic_algorithm) equation for LRS, I wrote the code mostly off the wikipedia article
// (1/n) * (sp - (2sp - 2)(i-1/n-1))
void linearRankSelection(vector<Trip>& genePool, vector<Trip>& parents,int selectionPressure) {
	double sp = double(selectionPressure);
	double n = static_cast<double>(genePool.size());	
	double firstFactor  = pow(n,-1); //n/1
	double secondFactor = 2 * sp - 2;

	//Mind you, we are creating a cumulative probability, so the random number will correspond to a range
	double p = firstFactor * sp;
	genePool[0].setLRSProb(p);
	for (int i = 1; i < n; i++) {
		double p = firstFactor * (sp - secondFactor * (((i + 1) - 1.0f) / (n - 1))) + genePool[i-1].getLRSProb();
		genePool[i].setLRSProb(p);
	}

	for (int i = 0; i < n/2; i++) {
		double a = genRandomDouble();
		Trip* previousGene = nullptr;
		for (auto& gene : genePool) {
			if (gene.getLRSProb() >= a) {
				if (previousGene != nullptr) {
					parents.push_back(*previousGene); //adding that gene to parents
					break;
				}
			}
			previousGene = &gene;
		}
	}
}

//selects the genes with the largest divergence from the mean of the current gene pool
void meanDevSelection(vector<Trip>& genePool, vector<Trip>& parents) {
	int n = static_cast<int>(genePool.size());
	float oneOverN = 1.0f / static_cast<double>(genePool.size());
	float fHat = 0;
	//calculate mean
	for (auto& gene : genePool) {
		fHat += gene.getPathLength();
	}
	float mean = oneOverN * fHat;
	
	for (auto& gene : genePool) {
		gene.setLRSProb(abs(gene.getPathLength() - mean));
	}
	//sort by largest divergence from mean
	sort(genePool.begin(), genePool.end(), compareMeanDevs);
	for (int i = 0; i < n / 2; i++) {
		parents.push_back(genePool[i]);
	}
}

//This maybe doesn't work - me, unenlightened
//It works now - me, so very enlightened
//Mind you this function does the same thing as RWS, kinda happy about this one
//Need to check how many gene's we are actually adding, set it to half pop size like
//I did for LRS
void newRWSSelection(vector<Trip>& genePool, vector<Trip>& parents, int popSize) {
	double n = static_cast<double>(genePool.size());
	float S = 0;

	//total sum of inverted paths
	for (auto& gene : genePool) {
		S = S + pow(gene.getPathLength(),-1);
	}

	for (int i = 0; i < n/2; i++) {
		float a = genRandom();
		float iSum = 0;
		int j = 0;
		do{
			iSum = iSum + (pow(genePool[j].getPathLength(),-1) / S);
			j++;
		} while (iSum < a && j < n-1);
		parents.push_back(genePool[j]);
	}

}

//original Roulette wheel, performs a bit better than the new one - no while loops
void RWS(vector<Trip>& genePool, vector<Trip>& parents, int popSize,float crossoverPer) {
	//setting up for the roulette wheel, need a total sum of inverted path's (inverted because we are minimizing not maximizing)
	float totalSumOfInvertedPaths = 0;
	for (auto& gene : genePool) {
		gene.calcInvertedProb(); //calculates the inverted path size for that path size
		totalSumOfInvertedPaths += gene.getInvertedPathLength();
	}

	//setting roulette prob - in other words, the inverted path length normalized to the sum of all the inverted paths
	for (auto& gene : genePool) {
		float rouletteProb = gene.getInvertedPathLength() / totalSumOfInvertedPaths; // = (1/pathLength) / totalSumOfInvertedPaths
		gene.setRouletteProb(rouletteProb);
	}

	//calculating the cumulative sum for each path in the gene pool
	float totalSum = 0;
	for (auto& gene : genePool) {
		totalSum += gene.getRouletteProb();
		gene.setCumProb(totalSum);
	}

	//creating a parents list
	for (int i = 0; i < int(crossoverPer * popSize); i++) { //iterating based on crossover per
		float rouletteRand = genRandom();
		for (auto& gene : genePool) {
			if (rouletteRand > gene.getCumProb()) { //if the random number is greater than that gene's cumulative prob, set to 1
				gene.setBoolProb(1);
			}
		}
		Trip* previousGene = nullptr;
		for (auto& gene : genePool) {
			if (gene.getBoolProb() == 0) { //finding last gene in cumulative prob less than random number
				if (previousGene != nullptr) {
					parents.push_back(*previousGene); //adding that gene to parents
					break;
				}
			}
			previousGene = &gene;
		}
		for (auto& gene : genePool) {
			gene.setBoolProb(0); //re setting back to 0
		}
	}
}


//picks a random index to splice genes, swaps parts, then changes anything that occurs twice
void spCrossover(Trip& gene1, Trip& gene2, vector<Trip>& children,int numCities) {
	int cut = rand() % (numCities - 1); //random placement 1-9 in this case
	vector<City> path1 = gene1.getPath(); //path 1
	vector<City> path2 = gene2.getPath(); //path 2
	vector<City> child1;
	vector<City> child2;

	//parent 1&2; 0:cut
	unordered_set<string> idsInChild2;
	unordered_set<string> idsInChild1;
	for (int i = 0; i <= cut; i++) {
		child1.push_back(path1[i]);
		child2.push_back(path2[i]);
		idsInChild1.insert(path1[i].getID());
		idsInChild2.insert(path2[i].getID());
	}

	for (auto& city : path1) {
		if (idsInChild2.insert(city.getID()).second) {  // Only inserts if ID is unique
			child2.push_back(city);
		}
	}

	for (auto& city : path2) {
		if (idsInChild1.insert(city.getID()).second) {  // Only inserts if ID is unique
			child1.push_back(city);
		}
	}

	Trip child1Trip(child1);
	Trip child2Trip(child2);
	child1Trip.calcPathLength(); //setting the path length
	child2Trip.calcPathLength();
	children.push_back(child1Trip);
	children.push_back(child2Trip);
}

//Uniform crossover works by giving a 20% for a swap between cities at the same index
//Once Swapped we prune over the entire gene to see which value we need to replace so there are no repeats
//Kind of slow if the chance is high, but works pretty well and is faster than spCrossover
void uniformCrossover(Trip& gene1, Trip& gene2, vector<Trip>& children,int numCities) {
	//slow but works, could find a way to improve the speed
	vector<City> child1 = gene1.getPath(); //path 1
	vector<City> child2 = gene2.getPath(); //path 2
	for (int i = 0; i < numCities; i++) {
		float chance = genRandom();
		City temp = child1[i];
		City temp2 = child2[i];
		if (chance > (1 - 0.2)) { //20%chance
			child1[i] = temp2;
			child2[i] = temp; //swapped
			int indexVal = i;
			for (int j = 0; j < numCities; j++) { //replace
				if (child1[j] == temp2 && j != indexVal) { //not index val makes sure we dont replace what we just swapped
					child1[j] = temp;
				}
				if (child2[j] == temp && j != indexVal) {
					child2[j] = temp2;
				}
			
			}
		}
	}

	Trip child1Trip(child1);
	Trip child2Trip(child2);
	child1Trip.calcPathLength(); //setting the path length
	child2Trip.calcPathLength();
	children.push_back(child1Trip);
	children.push_back(child2Trip);
}

//unlike the other crossover functions, this will only return one new unique gene, therefore, we adjust accordingly
//as far as I can tell this is working properly on the original dataset
//based on algorithm described here: https://en.wikipedia.org/wiki/Edge_recombination_operator
void edgeRecombination(Trip& gene1, Trip& gene2, vector<Trip>& children, int numCities) {
	vector<City> newGene;
	vector<City> path1 = gene1.getPath(); //path 1
	vector<City> path2 = gene2.getPath(); //path 2
	vector<vector<City>> g1Adj;
	vector<vector<City>> g2Adj;
	vector<vector<City>> unionList;
	map<City, vector<City>> mapping1;
	map<City, vector<City>> mapping2;

	//Create adjacency list
	g1Adj.push_back({ path1[numCities-1], path1[1] });
	g2Adj.push_back({ path2[numCities-1], path2[1] });
	for (int i = 1; i < numCities-1; i++) {
		g1Adj.push_back({ path1[i - 1], path1[i+1] });
		g2Adj.push_back({ path2[i - 1], path2[i+1] });
	}
	g1Adj.push_back({ path1[numCities - 2], path1[0] });
	g2Adj.push_back({ path2[numCities - 2], path2[0] });

	//mapping city to its adjacency list
	for (int i = 0; i < numCities; i++) {
		mapping1.insert({ path1[i], g1Adj[i] });
		mapping2.insert({ path2[i], g2Adj[i] });
	}

	//merge lists - all adjacent cities to same city in either path
	map<City, vector<City>> megaMap;
	for (int i = 0; i < numCities; i++) {
		City compCity = path1[i]; // Or path2[i], as they should represent the same cities.
		vector<City> combinedAdjList = mapping1[compCity];
		combinedAdjList.insert(
			combinedAdjList.end(),
			mapping2[compCity].begin(),
			mapping2[compCity].end()
		);
		//delete duplicates
		std::sort(combinedAdjList.begin(), combinedAdjList.end());
		combinedAdjList.erase(
			std::unique(combinedAdjList.begin(), combinedAdjList.end()),
			combinedAdjList.end()
		);

		megaMap[compCity] = combinedAdjList;
	}


	//Main Algorithm
	City N = path1[0];
	while (newGene.size() < numCities) {
		newGene.push_back(N);
		string id = N.getID();
		
		//remove N from ALL neighbor lists
		map<City, vector<City>>::iterator itt;
		for (itt = megaMap.begin(); itt != megaMap.end(); itt++) {
			City key = itt->first;
			vector<City> value = itt->second;
			value.erase(
				std::remove_if(
					value.begin(),
					value.end(),
					[id](City& city) {return city.getID() == id; }
				), value.end());
			megaMap[key] = value;
		}

		//Find neighbor of N with the smallest Number of neighbors itself
		//chooses a random one if more than one with the same size list
		int smallestSize = 999999;
		City nextN;
		vector<City> equalSmallest;
		if (megaMap[N].size() > 0) {
			for (auto city : megaMap[N]) {
				int adjSize = megaMap[city].size();
				if (adjSize < smallestSize) {
					equalSmallest.push_back(city);
				}
				else if (adjSize == smallestSize) {
					equalSmallest.push_back(city);
				}
			}
			int randIndex = rand() % equalSmallest.size();
			nextN = equalSmallest[randIndex];
		}
		else { //not used as much? i think
			//randomly chosen Node not in K = newGene, not super random, but much faster than the alternative I think
			for (const auto& city : path1) {
				// Check if city is not in newGene
				if (std::find(newGene.begin(), newGene.end(), city) == newGene.end()) {
					nextN = city;
					break; // Exit loop once we find a city not in newGene
				}
			}
		}

		N = nextN;
	}

	Trip childTrip(newGene);
	childTrip.calcPathLength(); //setting the path length
	children.push_back(childTrip);

}

//not able to implement this - but good I tried!
void PMX(Trip& gene1, Trip& gene2, vector<Trip>& children,int numCities) { //work in progress
	int left = rand() % (numCities - 3) + 1;
	int right_limit = numCities - 2;
	int middle = rand() % (right_limit - left) + left + 1;
	vector<City> path1 = gene1.getPath(); //path 1
	vector<City> path2 = gene2.getPath(); //path 2
	vector<City> child1;
	vector<City> child2;
	cout << "Left: " << left << " Middle: " << middle << endl;
	gene1.printPath(); cout << endl;
	gene2.printPath(); cout << endl;
	//populate children
	for (int i = 0; i < left; i++) { //left
		child1.push_back(path1[i]);
		child2.push_back(path2[i]);
	}
	for (int i = left; i <= middle; i++) { //middle - the key section
		child1.push_back(path2[i]);
		child2.push_back(path1[i]);
	}

	for (int i = middle + 1; i < numCities; i++) { //right
		child1.push_back(path1[i]);
		child2.push_back(path2[i]);
	}

	map<City, City> mapping1;
	map<City, City> mapping2;
	for (int i = left; i <= middle; i++) {
		mapping1.insert(pair<City, City>(child1[i], child2[i]));
		mapping2.insert(pair<City, City>(child2[i], child1[i]));
	}

	// Output the mappings
	map<City, City>::iterator ppt;
	cout << "Mapping 1:" << endl;
	for (ppt = mapping1.begin(); ppt != mapping1.end(); ppt++) {
		City key = ppt->first;
		City value = ppt->second;
		cout << key.getID() << ":" << value.getID() << "  ";
	}
	cout << endl;
	map<City, City>::iterator lpt;
	cout << "Mapping 2:" << endl;
	for (lpt = mapping2.begin(); lpt != mapping2.end(); lpt++) {
		City key = lpt->first;
		City value = lpt->second;
		cout << key.getID() << ":" << value.getID() << "  ";
	}
	cout << endl;

	map<City, City> mapping3;
	map<City, City> mapping4;

	map<City, City>::iterator it;
	map<City, City>::iterator jt;
	for (it = mapping1.begin(); it != mapping1.end(); it++) {
		City key = it->first;
		City value = it->second;
		for (jt = mapping1.begin(); jt != mapping1.end(); jt++) {
			City key1 = jt->first;
			City value1 = jt->second;
			if (key.getID() == value1.getID()) {
				mapping3.insert(pair<City, City>(value, key1));
			}
		}
	}

	map<City, City>::iterator ot;
	map<City, City>::iterator pt;
	for (ot = mapping2.begin(); ot != mapping2.end(); ot++) {
		City key = ot->first;
		City value = ot->second;
		for (pt = mapping2.begin(); pt != mapping2.end(); pt++) {
			City key1 = pt->first;
			City value1 = pt->second;
			if (key.getID() == value1.getID()) {
				mapping4.insert(pair<City, City>(value, key1));
			}
		}
	}

	map<City, City>::iterator zt;
	cout << "Mapping 3:" << endl;
	for (zt = mapping3.begin(); zt != mapping3.end(); zt++) {
		City key = zt->first;
		City value = zt->second;
		cout << key.getID() << ":" << value.getID() << "  ";
	}
	cout << endl;

	map<City, City>::iterator qt;
	cout << "Mapping 4:" << endl;
	for (qt = mapping4.begin(); qt != mapping4.end(); qt++) {
		City key = qt->first;
		City value = qt->second;
		cout << key.getID() << ":" << value.getID() << "  ";
	}
	cout << endl;



	cout << "end map" << endl;
	for (int i = 0; i < numCities; i++) {
		if (i > left-1 && i < middle + 1) {
			continue;
		}
		if (mapping3.find(child1[i]) != mapping3.end()) {
			child1[i] = mapping3[child1[i]];
		}
		else if (mapping1.find(child1[i]) != mapping1.end()) {
			child1[i] = mapping1[child1[i]];
		}

		if (mapping4.find(child2[i]) != mapping4.end()) {
			child2[i] = mapping4[child2[i]];
		}
		else if (mapping2.find(child2[i]) != mapping2.end()) {
			child2[i] = mapping2[child2[i]];
		}
	}

	cout << endl << "Child 1: " << endl;
	for (auto& city : child1) {
		cout << city.getID() << " ";
	}
	cout << endl;
	cout << endl;
	cout << "Child 2: " << endl;
	for (auto& city : child2) {
		cout << city.getID() << " ";
	}
	cout << endl;
	cout << endl;

}
