// openGL libraries
#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <GL/glut.h>

// armadillo library
#include <armadillo>

// c++ libraries
#include <map>


// my libraries
#include "Shape.h"

using namespace std;
using namespace arma;

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

void Shape::write_to_file()
{
  ofstream writer("topologycal_points_plane_medium.txt");
  int N = points.size();

  writer<<N<<endl;
  cout<<"we are about to begin"<<endl;

  mat S(N, N, fill::zeros);

  map<int, Vec3f> tag_points;
  map<Vec3f, int> points_tag;
  int counter = 0;
  for(auto itr = points.begin(); itr != points.end(); itr++)
  {
    //cout<<"we are writting!"<<endl;
    std::pair<Vec3f, int> c_1((*itr), counter);
    writer<<counter<<" "<<(*itr)<<endl;
    cout<<c_1.first<<endl;
    points_tag.insert(c_1);
    counter++;
  }

  writer<<segments.size()<<endl;
  for(auto itr = segments.begin(); itr != segments.end(); itr++)
  {
    cout<<"we are writting!"<<endl;
    int x = points_tag[(itr)->x];
    int y = points_tag[(itr)->y];
    writer<<x<<" "<<y<<endl;
  }

  writer<<faces.size()<<endl;
  for(auto itr = faces.begin(); itr != faces.end(); itr++)
  {
    cout<<"we are writting!"<<endl;
    int x = points_tag[(itr)->x];
    int y = points_tag[(itr)->y];
    int z = points_tag[(itr)->z];
    writer<<x<<" "<<y<<" "<<z<<endl;
  }

  writer.close();
}

void Shape::print()
{
  cout<<"Points:{ ";
  for(set<Vec3f>::iterator itr= points.begin(); itr != points.end(); itr++)
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

Shape& Shape::operator= (const Shape& obj)
{
  this->points = obj.points;
	this->segments = obj.segments;
	this->faces = obj.faces;
  return *this;
}

/*Shape Shape::operator+(const Shape& obj)
{
    Shape temp;
    temp.x = x + obj.x;
    temp.y = y + obj.y;
    temp.z = z + obj.z;
    return temp;
}

Shape Shape::translate(const Vec3f t)
{

}*/
