#include "../include/Particle.h"
//#include "../include/Vector.h"
//#include "../include/Pair.h"

//typedef Vector<double> Vec3D;
using namespace std;

Particle::Particle(vec p, vec vel, double mess, vec force)//: pos(p), v(vel), m(mess), F(force)
{
	this->pos = p;
	this->v = vel;

	this->F_list.resize(3);
	F_list[2] = vec(3, fill::zeros);
	F_list[2] = force;
	F_list[1] = vec(3, fill::zeros);
	vec gra(3, fill::zeros);
	gra(1) = -10 * mess;
	F_list[0] = vec(3, fill::zeros);
	F_list[0] = gra;
	fr_value = 0;

	this->pos_0 = this->pos;

	m = mess;

	construct_J();


}


void Particle::construct_J()
{
	J = Mat<double>(3, 12, fill::zeros);
	J(0,0) = 1;
	J(1,1) = 1;
	J(2,2) = 1;


	J(0,3) = pos(0);
	J(1,4) = pos(0);
	J(2,5) = pos(0);

	J(0,6) = pos(1);
	J(1,7) = pos(1);
	J(2,8) = pos(1);

	J(0,9) = pos(2);
	J(1,10) = pos(2);
	J(2,11) = pos(2);

	J_T = J.t();

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

bool operator< (const Particle& p1, const Particle& p2)
{
	if(p1.pos[0] < p2.pos[0])
		return true;
	if(p1.pos[0] == p2.pos[0] && p1.pos[1] < p2.pos[1])
		return true;
	if(p1.pos[0] == p2.pos[0] && p1.pos[1] == p2.pos[1] && p1.pos[2] < p2.pos[2])
		return true;
	return false;
}
