#include "City.h"

City::City()
{
	x = 0;
	y = 0;
	ID = "None";
}

City::City(float x, float y, string ID) : x(x),y(y),ID(ID)
{
}

float City::getX()
{
	return x;
}

float City::getY()
{
	return y;
}

string City::getID()
{
	return ID;
}

void City::setX(float newX)
{
	x = newX;
}

void City::setY(float newY)
{
	y = newY;
}

void City::setID(string newID)
{
	ID = newID;
}
