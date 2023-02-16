#include "../include/Face.h"

using namespace std;

Face::Face(Particle* x, Particle* y, Particle* z)
{
  this->x = x;
  this->y = y;
  this->z = z;
}

ostream& operator<<(ostream& os, const Face& obj)
{
  os<<"The x side of the face is:"<<endl<<*obj.x<<endl;
  os<<"The y side of the face is:"<<endl<<*obj.y<<endl;
  os<<"The z side of the face is:"<<endl<<*obj.z<<endl;
  return os;
}
