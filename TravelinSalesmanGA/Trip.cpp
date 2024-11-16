#include "Trip.h"
#include <math.h>

Trip::Trip()
{
	path;//init with an empty path
	pathSize = 0;
	fitnessProb = 0;
	cumProb = 0;
	rouletteProb = 0;
	boolProb = 0;
	lrsProb = 0;
	ersProb = 0;
}

Trip::Trip(vector<City> path): path(path)
{
	pathSize = 0;
	fitnessProb = 0;
	cumProb = 0;
	rouletteProb = 0;
	boolProb = 0;
	lrsProb = 0;
	ersProb = 0;
}

float Trip::getDistance(City city1, City city2)
{
	float x1 = city1.getX();
	float y1 = city1.getY();
	float x2 = city2.getX();
	float y2 = city2.getY();
	float distance = abs(sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2)));
	return distance;
}

void Trip::calcPathLength()
{
	pathSize = 0;
	for (auto i = path.begin();i != path.end(); ++i) {
		if (i + 1 != path.end()) {
			pathSize += getDistance(*i, *(i + 1));
		}
		else {
			pathSize += getDistance(*i, path.front());
		}
	}
}

void Trip::addCity(City city)
{
	path.push_back(city);
}

void Trip::swapCities(City city1, City city2)
{
}

vector<City> Trip::getPath()
{
	return path;
}

void Trip::printPath()
{
	for (auto& element : path) {
		std::cout << element.getID() << " ";
	}
	//std::cout << path.front().getID() << " ";
}

void Trip::calcInvertedProb()
{
	fitnessProb = 1 / (pathSize + 1e-10);
}

void Trip::setCumProb(float newProb)
{
	cumProb = newProb;
}

float Trip::getCumProb()
{
	return cumProb;
}

float Trip::getPathLength()
{
	return pathSize;
}

float Trip::getInvertedPathLength()
{
	return fitnessProb;
}

void Trip::setRouletteProb(float prob)
{
	rouletteProb = prob;
}

float Trip::getRouletteProb()
{
	return rouletteProb;
}


