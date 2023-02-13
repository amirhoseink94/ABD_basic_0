// openGL libraries
#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <GL/glut.h>

// my libraries
#include "../include/Shape.h"

using namespace std;

void Shape::draw_shape()
{
  for(set<Segment>::iterator itr=segments.begin(); itr!=segments.end(); itr++)
  {
    glBegin(GL_LINES);
    glVertex3f((*itr).x.x, (*itr).x.y, (*itr).x.z);
    glVertex3f((*itr).y.x, (*itr).y.y, (*itr).y.z);
    glEnd();
  }
}

void Shape::print()
{
  cout<<"Points:{ ";
  for(set<Vec3D>::iterator itr= points.begin(); itr != points.end(); itr++)
    cout<<(*itr);
  cout<<"}"<<endl;
  cout<<"Segments:{ ";
  for(set<Segment>::iterator itr= segments.begin(); itr != segments.end(); itr++)
    cout<<(*itr);
  cout<<"}"<<endl;
  cout<<"Faces:{ ";
  for(set<Face>::iterator itr= faces.begin(); itr != faces.end(); itr++)
    cout<<(*itr);
  cout<<"}"<<endl;
}
