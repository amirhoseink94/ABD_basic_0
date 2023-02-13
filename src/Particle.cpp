#include "../include/Particle.h"
#include "../include/Vector.h"
#include "../include/Pair.h"

//typedef Vector<float> Vec3D;

Particle::Particle(Vec3D p, Vec3D vel, float mess, Vec3D force)
{
	pos.x = p.x;
	pos.y = p.y;
	pos.z = p.z;
	
	v.x = vel.x;
	v.y = vel.y;
	v.z = vel.z;
	
	m = mess;
	F.x = force.x;
	F.y = force.y;
	F.z = force.z;
}



