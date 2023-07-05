#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <armadillo>
#include <iostream>

//#include "Body.h"
//#include "Segment.h"

using namespace std;
using namespace arma;
//using namespace ABD;


class Particle
{
	public:
	vec pos_0;
	vec v;
	double m;
	vector<vec> F_list;
	double fr_value;
	/*
		F[0] = graviti --fixed
		F[1] = friction
		F[2] = any outside force --so far we dont have any
	*/
	vec F;
	//Body* b;

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
