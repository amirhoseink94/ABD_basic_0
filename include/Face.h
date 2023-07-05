#ifndef FACE_H
#define FACE_H

// c++ libraries
#include <vector>
#include <armadillo>
#include <iostream>

// my libraries
#include "Particle.h"


// namespaces
using namespace std;
using namespace arma;

// typedef
typedef Mat<double> Vec3f;

class Face
{
  public:
  Particle* x;
  Particle* y;
  Particle* z;

  vec n;
  vec n_AB;
  vec n_BC;
  vec n_CA;


  double R;
  vec C_0;
  vec C;
  mat J_C;


  double r;
  vec c_0;
  vec c;
  mat J_c;

  Face(Particle* x, Particle* y, Particle* z);
  void Cal_centre();
  vector<Face*> nei_segments;
  friend ostream& operator<<(ostream& os, const Face& obj);
  bool operator==(const Face& s2) const;

  Face& operator= (const Face& obj)
	{
		this->x = obj.x;
		this->y = obj.y;
    this->z = obj.z;
		//this->m = obj.m;
		//this->F = obj.F;
    return *this;
	}
};

#endif
