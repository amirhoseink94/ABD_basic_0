#include "../include/Particle.h"
//#include "../include/Vector.h"
//#include "../include/Pair.h"

//typedef Vector<float> Vec3D;
using namespace std;

Particle::Particle(Vec3f p, Vec3f vel, float mess, Vec3f force)
{
	/*pos.x = p.x;
	pos.y = p.y;
	pos.z = p.z;

	v.x = vel.x;
	v.y = vel.y;
	v.z = vel.z;*/
	this->pos = p;
	this->v = vel;
	this->F = force;
	m = mess;
	/*F.x = force.x;
	F.y = force.y;
	F.z = force.z;*/
}

ostream& operator<<(ostream& os, const Particle& obj)
{
	os<<"================"<<endl;
	os<<"The particle position:"<<endl<<obj.pos<<endl;
	os<<"The particle velocity:"<<endl<<obj.v<<endl;
	os<<"The particle mess:"<<endl<<obj.m<<endl;
	os<<"The particle force:"<<endl<<obj.F<<endl;
	os<<"----------------"<<endl;
	return os;
}
