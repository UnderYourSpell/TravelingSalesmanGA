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

using namespace std;
using namespace std::chrono;

//Problem parameters
const int NUM_CITIES = 10;
const int POP_SIZE = 12;
const float CROSSOVER_PER = 0.5;
const float MUTATION_PER = 0.5; //50% mutation rate
const int ELITISM = 1;
const int REST = 10;
const int MAX_GENERATIONS = 100;

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

struct TSPProblemData {
	string name;
	string comment;
	string type;
	int dimension;
	string edge_weight_type;
	vector<City> initCities;
};

TSPProblemData readTSPFile(const string& filepath) {
	TSPProblemData data;
	ifstream infile(filepath);
	string line;

	if (!infile) {
		cerr << "Unable to open File" << filepath << endl;
		exit(1);
	}

	while (getline(infile, line)) {
		if (line.find("NAME") == 0) {
			data.name = line.substr(line.find(":") + 1);
		}
		else if (line.find("COMMENT") == 0) {
			data.comment = line.substr(line.find(":") + 1);
		}
		else if (line.find("TYPE") == 0) {
			data.type = line.substr(line.find(":") + 1);
		}
		else if (line.find("DIMENSION") == 0) {
			data.dimension = std::stoi(line.substr(line.find(":") + 1));
		}
		else if (line.find("EDGE_WEIGHT_TYPE") == 0) {
			data.edge_weight_type = line.substr(line.find(":") + 1);
		}
		else if (line == "NODE_COORD_SECTION") {
			string identifier;
			float x, y;
			while (infile >> identifier >> x >> y) {
				data.initCities.emplace_back(x, y, identifier);
			}
		}
	}
	return data;
}

int main() {
	srand(time(NULL));
	int run = 1;
	int roulette_wheel = 1; //use roulette wheel or not
	string filePath = "./tsp/original10.tsp";
	TSPProblemData data = readTSPFile(filePath);
	cout << data.name << endl;
	cout << data.comment << endl;
	cout << data.type << endl;
	cout << data.dimension << endl;
	cout << data.edge_weight_type << endl;
	vector<City> initCities = data.initCities;

	if (run == 0) {
		return 0;
	}
	else {
		auto start = high_resolution_clock::now();
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
			newTrip.calcPathLength();
			genePool.push_back(newTrip);
		}
	/*
	for (auto& gene : genePool) {
		gene.printPath();
		gene.printPathLength();
		cout << endl;
	}
	*/
		vector<Trip> newGen;
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
		}

		cout << endl << "Best Solution: " << endl;
		genePool[0].printPath();
		genePool[0].printPathLength();
		cout << endl;

		/*
		cout << "Final Generation" << endl;
		for (auto& gene : genePool) {
			cout << gene.getPathLength() << endl;
		}
		*/
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);
		cout << endl << "Time taken by function: "
			<< duration.count() << " micro Seconds" << endl;
		return 0;
	}
}