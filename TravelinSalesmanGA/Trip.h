#pragma once
#include <iostream>
#include "City.h"
#include <vector>
class Trip
{
public:
	Trip();
	Trip(vector<City> path);
	void addCity(City city);
	void swapCities(City city1, City city2);
	vector<City> getPath();
	void printPath();
	void calcPathLength();
	float getDistance(City city1, City city2);
	void printPathLength() {
		std::cout << pathSize;
	}
	void setPath(vector<City> newPath) {
		path = newPath;
	}
	void calcInvertedProb();
	void setCumProb(float newProb);
	float getCumProb();

	float getPathLength();
	float getInvertedPathLength();
	void setRouletteProb(float prob);
	float getRouletteProb();
	void setBoolProb(int newBool) {
		boolProb = newBool;
	}
	int getBoolProb(){
		return boolProb;
	}
private:
	vector<City> path;
	float pathSize;
	float fitnessProb;
	float cumProb;
	float rouletteProb;
	int boolProb;
};

