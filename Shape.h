#ifndef SHAPE_H
#define SHAPE_H

#include <vector>
#include "Point.h"
#include <utility>
#include <tuple>

using namespace std;

class Shape
{
	public:
	vector<Point> points;
	vector< pair<Point, Point> > segments;
	vector< tuple<Point, Point, Point> > faces;	 
};

#endif
