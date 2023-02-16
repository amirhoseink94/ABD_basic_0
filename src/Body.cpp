#include "../include/Body.h"

#include <vector>
using namespace std;
using namespace arma;

void Body::update()
{
  for(long unsigned int i=0; i<points.size(); i++)
  {
    for(long unsigned int j=i+1; j<points.size(); j++)
    {
      if(approx_equal(points[i].pos, points[j].pos, "absdiff", 0.001)) //(points[i] == points[j])
      {
        points.erase(points.begin() + j);
        j--;
      }
    }
  }

  for(long unsigned int i=0; i<segments.size(); i++)
  {
    for(long unsigned int j=i+1; j<segments.size(); j++)
    {
      if(segments[i] == segments[j])
      {
        //remove thing
      }
    }
  }

  /*for(int i=0; i<faces.size(); i++)
  {
    for(int j=i+1; j<faces.size(); j++)
    {
      if(*faces[i] == *faces[j])
      {
        //remove thing
      }
    }
  }*/

}
