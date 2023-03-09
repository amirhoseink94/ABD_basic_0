#ifndef SHAPE_H
#define SHAPE_H

#include <set>
#include <iostream>
#include <stdio.h>
#include <bits/stdc++.h>

#include "Vector.h"
#include "Pair.h"

typedef Pair<Vec3f> Segment;
typedef Vector<Vec3f> Face;

using namespace std;

class Shape
{
	public:
	set<Vec3f> points;
	set<Segment> segments;
	set<Face> faces;

	void draw_shape();
	void print();
	Shape& operator= (const Shape& obj);

	void write_to_file();
	//Shape operator+(const Shape& obj);
};

#endif
