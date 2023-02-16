#include "../include/Utilities.h"

#include <armadillo>

using namespace std;
using namespace arma;

typedef Mat<float> Vec3f;

Body Util::make_mesh(Body start, int n)
{
    Body res = Util::make_mesh_handler(start, 0, n);
    return res;
}

Body Util::make_mesh_handler(Body start, int itr, int n)
{
  if(itr == n)
  {
    return start;
  }

  Vec3f Z(3,1, fill::zeros);
  Particle A(start.points[0].pos, Z, 1, Z);
  Particle B(start.points[1].pos, Z, 1, Z);
  Particle C(start.points[2].pos, Z, 1, Z);
  Particle m_AB((A.pos + B.pos)*(0.5), Z, 1, Z);
  Particle m_BC((B.pos + C.pos)*(0.5), Z, 1, Z);
  Particle m_CA((C.pos + A.pos)*(0.5), Z, 1, Z);
  //cout<<C<<A<<C + A<<endl;
  Segment m_AB_A(&m_AB, &A);
  Segment m_AB_B(&m_AB, &B);
  Segment m_BC_B(&m_BC, &B);
  Segment m_BC_C(&m_BC, &C);
  Segment m_CA_C(&m_CA, &C);
  Segment m_CA_A(&m_CA, &A);
  Segment m_AB_m_CA(&m_AB, &m_CA);
  Segment m_AB_m_BC(&m_AB, &m_BC);
  Segment m_BC_m_CA(&m_BC, &m_CA);

  Face m_AB_m_CA_A(&m_AB, &m_CA, &A);
  Face m_AB_m_BC_m_CA(&m_AB, &m_CA, &m_BC);
  Face m_AB_m_BC_B(&m_AB, &m_BC, &B);
  Face m_BC_m_CA_C(&m_BC, &m_CA, &C);

  Body up;
  up.points.push_back(A);
  up.points.push_back(m_AB);
  up.points.push_back(m_CA);

  up.segments.push_back(m_AB_A);
  up.segments.push_back(m_CA_A);
  up.segments.push_back(m_AB_m_CA);

  up.faces.push_back(m_AB_m_CA_A);
  // ----------------------------

  Body left;
  left.points.push_back(B);
  left.points.push_back(m_AB);
  left.points.push_back(m_BC);

  left.segments.push_back(m_AB_B);
  left.segments.push_back(m_BC_B);
  left.segments.push_back(m_AB_m_CA);

  left.faces.push_back(m_AB_m_BC_B);
  // ----------------------------

  Body right;
  right.points.push_back(C);
  right.points.push_back(m_BC);
  right.points.push_back(m_CA);

  right.segments.push_back(m_BC_C);
  right.segments.push_back(m_CA_C);
  right.segments.push_back(m_BC_m_CA);

  right.faces.push_back(m_BC_m_CA_C);
  // ----------------------------

  Body middle;
  middle.points.push_back(m_AB);
  middle.points.push_back(m_BC);
  middle.points.push_back(m_CA);

  middle.segments.push_back(m_AB_m_CA);
  middle.segments.push_back(m_BC_m_CA);
  middle.segments.push_back(m_AB_m_BC);

  middle.faces.push_back(m_AB_m_BC_m_CA);
  //----
  Body res_up = make_mesh_handler(up, itr+1, n);
  Body res_middle = make_mesh_handler(middle, itr+1, n);
  Body res_right = make_mesh_handler(right, itr+1, n);
  Body res_left = make_mesh_handler(left, itr+1, n);
  Body res = merge_shapes(merge_shapes(res_left, res_right), merge_shapes(res_up, res_middle));


  return res;
}

Body Util::merge_shapes(Body A, Body B)
{
  Body m;
  m.points = A.points;
  m.segments = A.segments;
  m.faces = A.faces;

  m.points.insert(m.points.end(), B.points.begin(), B.points.end());
  m.segments.insert(m.segments.end(), B.segments.begin(), B.segments.end());
  m.faces.insert(m.faces.end(), B.faces.begin(), B.faces.end());

  return m;
  /*int n_points = A.points.size();
  int n_segment = A.segments.size();
  int n_face = A.faces.size();

  for(int i=0; i<B.points.size(); i++)
  {
    bool check_find = true;
    for(int j=0; j<n_points; j++)
    {
      if(A.points[j] == B.points[i])
      {
        check_find =false;
        break;
      }
    }
    if(check_find)
    {
        m.points.push_back(B.points[i]);
    }
  }
  for(int i=0; i<B.segments.size(); i++)
  {
    bool check_find = true;
    for(int j=0; j<n_segments; j++)
    {
      if(A.segments[j] == B.segments[i])
      {
        check_find =false;
        break;
      }
    }
    if(check_find)
    {
        m.segments.push_back(B.segments[i]);
    }
  }
  for(int i=0; i<B.faces.size(); i++)
  {
    bool check_find = true;
    for(int j=0; j<n_faces; j++)
    {
      if(A.faces[j] == B.faces[i])
      {
        check_find =false;
        break;
      }
    }
    if(check_find)
    {
        m.faces.push_back(B.faces[i]);
    }
  }


  return m;*/
}


/*Shape Util::make_mesh(Shape start, int n)
{
    Shape res = Util::make_mesh_handler(start, 0, n);
    return res;
}

Shape Util::make_mesh_handler(Shape start, int itr, int n)
{
  if(itr == n)
  {
    return start;
  }
  //cout<<itr<<endl;
  set<Vec3D>::iterator set_itr = start.points.begin();
  Vec3D A((*set_itr).x, (*set_itr).y, (*set_itr).z);
  set_itr++;
  Vec3D B((*set_itr).x, (*set_itr).y, (*set_itr).z);
  set_itr++;
  Vec3D C((*set_itr).x, (*set_itr).y, (*set_itr).z);

  Vec3D m_AB = (A + B).num_multi(0.5);
  Vec3D m_BC = (B + C).num_multi(0.5);
  Vec3D m_CA = (C + A).num_multi(0.5);
  //cout<<C<<A<<C + A<<endl;
  Segment m_AB_A(m_AB, A);
  Segment m_AB_B(m_AB, B);
  Segment m_BC_B(m_BC, B);
  Segment m_BC_C(m_BC, C);
  Segment m_CA_C(m_CA, C);
  Segment m_CA_A(m_CA, A);
  Segment m_AB_m_CA(m_AB, m_CA);
  Segment m_AB_m_BC(m_AB, m_BC);
  Segment m_BC_m_CA(m_BC, m_CA);

  Face m_AB_m_CA_A(m_AB, m_CA, A);
  Face m_AB_m_BC_m_CA(m_AB, m_CA, m_BC);
  Face m_AB_m_BC_B(m_AB, m_BC, B);
  Face m_BC_m_CA_C(m_BC, m_CA, C);

  Shape up;
  up.points.insert(A);
  up.points.insert(m_AB);
  up.points.insert(m_CA);

  up.segments.insert(m_AB_A);
  up.segments.insert(m_CA_A);
  up.segments.insert(m_AB_m_CA);

  up.faces.insert(m_AB_m_CA_A);
  // ----------------------------

  Shape left;
  left.points.insert(B);
  left.points.insert(m_AB);
  left.points.insert(m_BC);

  left.segments.insert(m_AB_B);
  left.segments.insert(m_BC_B);
  left.segments.insert(m_AB_m_CA);

  left.faces.insert(m_AB_m_BC_B);
  // ----------------------------

  Shape right;
  right.points.insert(C);
  right.points.insert(m_BC);
  right.points.insert(m_CA);

  right.segments.insert(m_BC_C);
  right.segments.insert(m_CA_C);
  right.segments.insert(m_BC_m_CA);

  right.faces.insert(m_BC_m_CA_C);
  // ----------------------------

  Shape middle;
  middle.points.insert(m_AB);
  middle.points.insert(m_BC);
  middle.points.insert(m_CA);

  middle.segments.insert(m_AB_m_CA);
  middle.segments.insert(m_BC_m_CA);
  middle.segments.insert(m_AB_m_BC);

  middle.faces.insert(m_AB_m_BC_m_CA);
  //----
  Shape res_up = make_mesh_handler(up, itr+1, n);
  Shape res_middle = make_mesh_handler(middle, itr+1, n);
  Shape res_right = make_mesh_handler(right, itr+1, n);
  Shape res_left = make_mesh_handler(left, itr+1, n);
  Shape res = merge_shapes(merge_shapes(res_left, res_right), merge_shapes(res_up, res_middle));


  return res;
}

Shape Util::merge_shapes(Shape s, Shape v)
{
  Shape m;
  m.points = s.points;
  m.segments = s.segments;
  m.faces = s.faces;
  for(set<Vec3D>::iterator itr=v.points.begin(); itr != v.points.end(); itr++)
    m.points.insert((*itr));

  for(set<Segment>::iterator itr=v.segments.begin(); itr != v.segments.end(); itr++)
    m.segments.insert((*itr));

  for(set<Face>::iterator itr=v.faces.begin(); itr !=v.faces.end(); itr++)
    m.faces.insert((*itr));

  return m;
}*/
