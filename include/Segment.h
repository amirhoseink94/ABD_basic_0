#ifndef SEGMENT_H
#define SEGMENT_H

// c++ libraries
#include <vector>

// my libraries
#include "Particle.h"
#include "Face.h"

// namespasec
using namespace std;
using namespace arma;

//typedef


class Segment
{
  public:
  Particle* x;
  Particle* y;



	vector<Particle*> nei_particles;
	vector<Face*> nei_faces;
  Segment(Particle* x, Particle* y)
  {
    this->x = x;
    this->y = y;
  }
  bool operator==(const Segment& rhs) const;
  void print();
  friend ostream& operator<<(ostream& os, const Segment& obj);

  Segment& operator= (const Segment& obj)
	{
		this->x = obj.x;
		this->y = obj.y;
		//this->m = obj.m;
		//this->F = obj.F;

    return *this;
	}

};

#endif
