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
const int NUM_CITIES = 48;

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
