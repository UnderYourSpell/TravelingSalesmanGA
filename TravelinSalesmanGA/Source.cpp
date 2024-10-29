#include <iostream>
#include <string>
#include <vector>
#include "City.h"
#include "Trip.h"
#include <algorithm>
#include <chrono>
#include <numeric>
#include <time.h>

using namespace std;
using namespace std::chrono;

//Problem parameters
const int NUM_CITIES = 10;
const int POP_SIZE = 12;
const float CROSSOVER_PER = 0.5;
const float MUTATION_PER = 0.5; //20% mutation rate
const int ELITISM = 1;
const int REST = 10;
const int MAX_GENERATIONS = 500;

float genRandom() { //generates random number between 0 and 1
	return ((float)rand()) / RAND_MAX;
}

bool comparePaths(Trip i1, Trip i2) {
	return(i1.getPathLength() < i2.getPathLength());
}


void spCrossover(Trip& gene1, Trip& gene2, vector<Trip> &children) {
	int cut = rand() % (NUM_CITIES-1); //random placement 1-9 in this case
	vector<City> path1 = gene1.getPath(); //path 1
	vector<City> path2 = gene2.getPath(); //path 2
	vector<City> child1;
	vector<City> child2;
	
	//parent 1&2; 0:cut
	for (int i = 0; i <= cut; i++) {
		child1.push_back(path1[i]);
		child2.push_back(path2[i]);
	}
	//can consider hashing as well because this is slow

	//setting child 1
	for (size_t i = 0; i < path2.size(); i++) {
		int found = 0;
		for (size_t j = 0; j < child1.size(); j++) {
			if (child1[j].getID() == path2[i].getID()) {
				found = 1;
				break;
			}
		}
		if (found != 1) child1.push_back(path2[i]);
	}
	//setting child 2
	for (size_t i = 0; i < path1.size(); i++) {
		int found = 0;
		for (size_t j = 0; j < child2.size(); j++) {
			if (child2[j].getID() == path1[i].getID()) {
				found = 1;
				break;
			}
		}
		if (found != 1) child2.push_back(path1[i]);
	}
	Trip child1Trip(child1);
	Trip child2Trip(child2);
	child1Trip.calcPathLength(); //setting the path length
	child2Trip.calcPathLength();
	children.push_back(child1Trip);
	children.push_back(child2Trip);
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

int main() {
	srand(time(NULL));
	auto start = high_resolution_clock::now();

	float x[NUM_CITIES] = { 0, 3, 6, 7, 15, 10, 16, 5, 8, 1.5 };
	float y[NUM_CITIES] = { 1,2,1,4.5,-1,2.5,11,6,9,12 };
	string city_names[NUM_CITIES] = { "Gliwice", "Cairo", "Rome", "Krakow", "Paris", "Alexandria", "Berlin", "Tokyo", "Hong Kong", "Rio" };

	vector<City> initCities;
	for (int i = 0; i < NUM_CITIES; i++) {
		City newCity = City(x[i], y[i], city_names[i]);
		initCities.push_back(newCity);
	}

	vector<Trip> genePool;
	int picks[NUM_CITIES] = {};
	for (int i = 0; i < NUM_CITIES; i++) picks[i] = i; //populating picks

	for (int i = 0; i < POP_SIZE; i++) {
		Trip newTrip;
		int newPicks[NUM_CITIES];
		copy(begin(picks), end(picks), begin(newPicks));
		int n = static_cast<int>(sizeof(newPicks) / sizeof(*newPicks));
		//my way of populating random genes O(n)
		for (int i = 0; i < NUM_CITIES; i++) {
			int randIndex = rand() % n;
			int numToAdd = newPicks[randIndex];
			newPicks[randIndex] = newPicks[n - 1];
			n--;
			newTrip.addCity(initCities[numToAdd]);
		}
		genePool.push_back(newTrip);
	}

	for (auto& gene : genePool) {
		gene.calcPathLength();
		gene.printPath();
		gene.printPathLength();
		cout << endl;
	}

	vector<Trip> bestOfGen;

	for (int p = 0; p <= MAX_GENERATIONS; p++) {
		
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
		vector<Trip> parents;
		for (int i = 0; i < int(CROSSOVER_PER * POP_SIZE); i++){ //iterating based on crossover per
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

		//__ __ __ two per, iterating
		vector<Trip> children;
		for(size_t i = 0;i+1<parents.size();i+=2){
			spCrossover(parents[i], parents[i + 1],children);
		}
		//mutation  - swapping cities in a path - 20% mutation chance per gene in children pool - introducing new genes essentially
		for (size_t i = 0; i < children.size(); i++) {
			float mutateThreshold = genRandom();
			if (mutateThreshold > (1 - MUTATION_PER)) {
				mutate(children[i]);
			}
		}
		
		//now we add the children to the current population - sort - then move on
		//1. Grab top 2 from original gene pool
		//2. Put all parents and children in the same vector - add children to parent vector
		//3. sort the parent vector, grab enough to fill enough for a population
		vector<Trip> newGen;
		sort(genePool.begin(), genePool.end(), comparePaths); //sort the original
		//1.
		for (auto i = genePool.begin(); i != genePool.begin() + ELITISM; ++i) {
			newGen.push_back(*i);
		}
		//2.
		parents.insert(parents.end(), children.begin(), children.end());
		sort(parents.begin(), parents.end(), comparePaths); //sort the original
		//3.
		for (auto i = parents.begin(); i != parents.end()-ELITISM; ++i) {
			newGen.push_back(*i);
		}

		genePool.clear();
		genePool = newGen;
		newGen.clear();
		bestOfGen.push_back(genePool[0]);
	}

	bestOfGen.back().printPath();
	bestOfGen.back().printPathLength();
	

	/*
	cout << "Final Generation" << endl;
	for (auto& gene : genePool) {
		cout << gene.getPathLength() << endl;
	}
	*/
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout << endl << "Time taken by function: "
		<< duration.count() << " microseconds" << endl;
	return 0;
}