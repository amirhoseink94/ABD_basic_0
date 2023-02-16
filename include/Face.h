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
typedef Mat<float> Vec3f;

class Face
{
  public:
  Particle* x;
  Particle* y;
  Particle* z;

	vector<Face*> nei_segments;
  Face(Particle* x, Particle* y, Particle* z);
  friend ostream& operator<<(ostream& os, const Face& obj);
};

#endif
