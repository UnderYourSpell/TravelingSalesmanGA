#pragma once
#include <math.h>
#include "City.h"
#include "Trip.h"

inline float getDistance(City city1, City city2);

//creates a Trip based on the greedy nearest neighbors
//inline  - in header file
inline Trip NearestNeighbor(vector<City> initCities, int numCities) {
	Trip newTrip;
	vector<City> U; //new gene we will add to
	vector<City> V; //random path we build to start the process

	vector<int> picks;
	for (int i = 0; i < numCities; i++) { picks.push_back(i); };
	vector<int> newPicks;
	newPicks.insert(newPicks.begin(), picks.begin(), picks.end());
	int n = numCities;
	//create a new random path
	for (int i = 0; i < numCities; i++) {
		int randIndex = rand() % n;
		int numToAdd = newPicks[randIndex];
		newPicks[randIndex] = newPicks[n - 1];
		n--;
		V.push_back(initCities[numToAdd]);
	}

	//Using U and V, U is the nearest neighbor path, V is the original path
	U.push_back(V[0]);
	V.erase(V.begin());

	//While U has less than the number of cities, add the nearest city to the last city added to V
	do {
		float dist = 999999999999999; 
		City closestCity;
		int deleteIndex = 0;
		//find min distance from U[0] to whatever city in v
		for (int i = 0; i < V.size(); i++) {
			float curDist = getDistance(U.back(), V[i]);
			if (curDist < dist) {
				closestCity = V[i];
				dist = curDist;
				deleteIndex = i;
			}
		}
		U.push_back(closestCity);
		V.erase(V.begin() + deleteIndex);
	} while (U.size() < numCities);

	for (auto& city : U) {
		newTrip.addCity(city);
	}
	newTrip.calcPathLength();
	return newTrip;
}

inline float getDistance(City city1, City city2)
{
	float x1 = city1.getX();
	float y1 = city1.getY();
	float x2 = city2.getX();
	float y2 = city2.getY();
	float distance = sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2));
	return distance;
}
