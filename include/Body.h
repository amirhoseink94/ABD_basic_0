#ifndef BODY_H
#define BODY_H

#include <set>
#include <iostream>
#include <stdio.h>
#include <armadillo>
#include<bits/stdc++.h>

#include "Particle.h"
#include "Segment.h"
#include "Face.h"


using namespace std;
using namespace arma;

typedef Mat<float> Vec3f;

class Body
{
	public:
	vector<Particle> points;
	vector<Segment> segments;
	vector<Face> faces;

	void update();

	//void draw_body();
	//void print();

};

#endif
