#include "../include/Segment.h"
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;
bool Segment::operator==(const Segment& s2) const
{
  if(approx_equal(x->pos, s2.x->pos, "absdiff", 0.001) && approx_equal(y->pos, s2.y->pos, "absdiff", 0.001))
    return true;
  if(approx_equal(x->pos, s2.y->pos, "absdiff", 0.001) && approx_equal(y->pos, s2.x->pos, "absdiff", 0.001))
    return true;
  return false;
}

ostream& operator<<(ostream& os, const Segment& obj)
{
  os<<"The x side of the segment is:"<<endl<<*obj.x<<endl;
  os<<"The y side of the segment is:"<<endl<<*obj.y<<endl;
  return os;
}
