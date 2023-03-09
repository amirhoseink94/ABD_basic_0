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


bool Face::operator==(const Face& s2) const
{
  if(approx_equal(x->pos, s2.x->pos, "absdiff", 0.001) && approx_equal(y->pos, s2.y->pos, "absdiff", 0.001) && approx_equal(z->pos, s2.z->pos, "absdiff", 0.001))
    return true;
  if(approx_equal(x->pos, s2.x->pos, "absdiff", 0.001) && approx_equal(y->pos, s2.z->pos, "absdiff", 0.001) && approx_equal(z->pos, s2.y->pos, "absdiff", 0.001))
    return true;
  if(approx_equal(x->pos, s2.y->pos, "absdiff", 0.001) && approx_equal(y->pos, s2.x->pos, "absdiff", 0.001) && approx_equal(z->pos, s2.z->pos, "absdiff", 0.001))
    return true;
  if(approx_equal(x->pos, s2.y->pos, "absdiff", 0.001) && approx_equal(y->pos, s2.z->pos, "absdiff", 0.001) && approx_equal(z->pos, s2.x->pos, "absdiff", 0.001))
    return true;
  if(approx_equal(x->pos, s2.z->pos, "absdiff", 0.001) && approx_equal(y->pos, s2.x->pos, "absdiff", 0.001) && approx_equal(z->pos, s2.y->pos, "absdiff", 0.001))
      return true;
  if(approx_equal(x->pos, s2.z->pos, "absdiff", 0.001) && approx_equal(y->pos, s2.y->pos, "absdiff", 0.001) && approx_equal(z->pos, s2.x->pos, "absdiff", 0.001))
    return true;
  return false;
}
