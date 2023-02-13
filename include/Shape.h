#ifndef SHAPE_H
#define SHAPE_H

#include <set>
#include <iostream>
#include <stdio.h>
#include<bits/stdc++.h>

#include "Vector.h"
#include "Pair.h"

typedef Pair<Vec3D> Segment;
typedef Vector<Vec3D> Face;

using namespace std;

class Shape
{
	public:
	set<Vec3D> points;
	set<Segment> segments;
	set<Face> faces;

	void draw_shape();
	void print();
};

#endif
