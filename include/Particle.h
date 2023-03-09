#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <armadillo>
#include <iostream>

//#include "Face.h"
//#include "Segment.h"

using namespace std;
using namespace arma;


class Particle
{
	public:
	vec pos_0;
	vec v;
	double m;
	vec F;

	vec pos;

	Mat<double> J;
	Mat<double> J_T;

//	vector<Segment*> nei_segments;
//	vector<Face*> nei_faces;

	Particle(vec, vec, double, vec); //pos - vel - mess - force


	Particle& operator= (const Particle& obj)
	{
		this->pos = obj.pos;
		this->v = obj.v;
		this->m = obj.m;
		this->F = obj.F;

		return *this;
	}
	friend ostream& operator<<(ostream& os, const Particle& obj);
	friend bool operator< (const Particle& p1, const Particle& p2);

	private:
	void construct_J();
};

#endif
