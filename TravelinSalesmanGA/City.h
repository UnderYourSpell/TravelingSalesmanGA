#pragma once
#include <string>
using namespace std;

class City
{
public:
	City();
	City(float x,float y,string ID);
	float getX();
	float getY();
	string getID();
	void setX(float newX);
	void setY(float newY);
	void setID(string ID);
	bool operator<(const City& other) const {
		return ID < other.ID;
	}
	bool operator!=(const City& other) const {
		return ID != other.ID;
	}
	bool operator==(const City& other) const {
		return ID == other.ID;
	}
private:
	float x;
	float y;
	string ID;
};

