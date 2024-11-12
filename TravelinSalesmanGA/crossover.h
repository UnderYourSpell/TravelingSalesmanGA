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

using namespace std;
const int NUM_CITIES = 10;

float genRandom() { //generates random number between 0 and 1
	return ((float)rand()) / RAND_MAX;
}

void mutate(Trip& gene) { //swapping two random cities
	int firstSwapIndex = rand() % (NUM_CITIES - 1);
	int secondSwapIndex = rand() % (NUM_CITIES - 1);
	vector<City> mutatedPath = gene.getPath();
	City temp = mutatedPath[firstSwapIndex];
	mutatedPath[firstSwapIndex] = mutatedPath[secondSwapIndex];
	mutatedPath[secondSwapIndex] = temp;
	gene.setPath(mutatedPath);
	gene.calcPathLength();
}

void scrambleMutate(Trip& gene,int mutateLength) {
	vector<City> path = gene.getPath();
	int startIndex = rand() % (NUM_CITIES - mutateLength - 1);
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

void newRWSSelection(vector<Trip>& genePool, vector<Trip>& parents, int popSize) {
	int n = genePool.size();
	float S = 0;

	for (auto& gene : genePool) {
		S += 1.0f / gene.getPathLength();
	}
	cout << "S:" << S;
	for (int i = 0; i < popSize; i++) {
		float a = static_cast<float>(rand()) / RAND_MAX * S;
		float iSum = 0;
		int j = 0;
		do{
			iSum = iSum + (1.0f / genePool[j].getPathLength());
			j++;
		} while (iSum < a && j < n-1);
		parents.push_back(genePool[j]);
	}
}

void RWS(vector<Trip>& genePool, vector<Trip>& parents, int popSize,float crossoverPer) {
	//setting up for the roulette wheel, need a total sum of inverted path's (inverted because we are minimizing not maximizing)
	float totalSumOfInvertedPaths = 0;
	for (auto& gene : genePool) {
		gene.calcInvertedProb(); //calculates the inverted path size for that path size
		totalSumOfInvertedPaths += gene.getInvertedPathLength();
	}

	//setting roulette prob - in other words, the inverted path length normalized to the sum of all the inverted paths
	for (auto& gene : genePool) {
		float rouletteProb = gene.getInvertedPathLength() / totalSumOfInvertedPaths;
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
void spCrossover(Trip& gene1, Trip& gene2, vector<Trip>& children) {
	int cut = rand() % (NUM_CITIES - 1); //random placement 1-9 in this case
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
void uniformCrossover(Trip& gene1, Trip& gene2, vector<Trip>& children) {
	//slow but works, could find a way to improve the speed
	vector<City> child1 = gene1.getPath(); //path 1
	vector<City> child2 = gene2.getPath(); //path 2
	for (int i = 0; i < NUM_CITIES; i++) {
		float chance = genRandom();
		City temp = child1[i];
		City temp2 = child2[i];
		if (chance > (1 - 0.2)) { //20%chance
			child1[i] = temp2;
			child2[i] = temp; //swapped
			int indexVal = i;
			for (int j = 0; j < NUM_CITIES; j++) {
				if (child1[j] == temp2 && j != indexVal) {
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




void partiallyMappedCrossover(Trip& gene1, Trip& gene2, vector<Trip>& children) { //work in progress
	int left = rand() % (NUM_CITIES - 3) + 1;
	int right_limit = NUM_CITIES - 2;
	int middle = rand() % (right_limit - left) + left + 1;
	vector<City> path1 = gene1.getPath(); //path 1
	vector<City> path2 = gene2.getPath(); //path 2
	vector<City> child1;
	vector<City> child2;
	cout << "Left: " << left << " Middle: " << middle << endl;


	//populate children
	for (int i = 0; i <= left; i++) { //left
		child1.push_back(path1[i]);
		child2.push_back(path2[i]);
	}

	for (int i = left + 1; i <= middle; i++) { //middle
		child1.push_back(path2[i]);
		child2.push_back(path1[i]);
	}

	for (int i = middle + 1; i < NUM_CITIES; i++) { //right
		child1.push_back(path1[i]);
		child2.push_back(path2[i]);
	}


	map<City, City> mapping1;
	map<City, City> mapping2;
	for (int i = left; i <= middle; i++) {
		mapping1.insert(pair<City, City>(child1[i], child2[i]));
		mapping2.insert(pair<City, City>(child2[i], child1[i]));
	}

	cout << "Mapping 1" << endl;
	for (int i = left; i <= middle; i++) {
		cout << child2[i].getID() << ":" << mapping1[child2[i]].getID() << endl;
	}
	cout << "Mapping 2" << endl;
	for (int i = left; i <= middle; i++) {
		cout << child1[i].getID() << ":" << mapping2[child1[i]].getID() << endl;
	}

	for (int i = 0; i < NUM_CITIES; i++) {
		// Skip the middle segment, which is already mapped
		if (i < left || i > middle) {
			// Legalize child1[i] with mapping1
			if (mapping1.count(child1[i]) > 0) {
				if (mapping1[child1[i]] != child1[i]) {
					while (mapping1.count(child1[i]) > 0) {
						cout << "skdjskdj" << endl;
						child1[i] = mapping1[child1[i]];
					}
				}
			}

			// Legalize child2[i] with mapping2
			if (mapping1.count(child2[i]) > 0) {
				if (mapping2[child2[i]] != child2[i]) {
					while (mapping1.count(child2[i]) > 0) {
						child2[i] = mapping2[child2[i]];
					}
				}
			}
		}
	}

	Trip child1Trip(child1);
	Trip child2Trip(child2);
	//child1Trip.calcPathLength(); //setting the path length
	//child1Trip.printPath();
	//child2Trip.printPath();
	cout << "Parent 1: ";
	gene1.printPath();
	Trip child1Trip1(child1);
	cout << "Child1    ";
	child1Trip1.printPath();
	cout << "Parent 2: ";
	gene2.printPath();
	cout << "Child2    ";
	Trip child1Trip2(child2);
	child1Trip2.printPath();
	child2Trip.calcPathLength();
	children.push_back(child1Trip);
	children.push_back(child2Trip);

}
