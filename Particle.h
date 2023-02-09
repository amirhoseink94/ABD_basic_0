#ifndef PARTICLE_H
#define PARTICLE_H

#include "Vector.h"

using namespace std;
typedef Vector<float> Vec3D;

class Particle
{	
	typedef Vector<float> Vec3D;
	public:
	Vec3D pos;
	Vec3D v;
	float m;
	Vec3D F;
	
	Particle(Vec3D, Vec3D, float, Vec3D);
};

#endif
