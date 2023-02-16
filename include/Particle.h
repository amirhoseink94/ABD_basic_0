#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <armadillo>
#include <iostream>

//#include "Face.h"
//#include "Segment.h"

using namespace std;
using namespace arma;
typedef Mat<float> Vec3f;

class Particle
{
	public:
	Vec3f pos;
	Vec3f v;
	float m;
	Vec3f F;

//	vector<Segment*> nei_segments;
//	vector<Face*> nei_faces;

	Particle(Vec3f, Vec3f, float, Vec3f); //pos - vel - mess - force

	friend ostream& operator<<(ostream& os, const Particle& obj);
};

#endif
