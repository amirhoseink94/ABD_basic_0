#include "../include/Utilities.h"

using namespace std;


Shape Util::make_mesh(Shape start, int n)
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
  /*cout<<"we are here!"<<endl;
  set<Vec3D>::iterator set_itr = start.points.begin();
  Vec3D A((*set_itr).x, (*set_itr).y, (*set_itr).z);
  set_itr++;
  Vec3D B((*set_itr).x, (*set_itr).y, (*set_itr).z);
  set_itr++;
  Vec3D C((*set_itr).x, (*set_itr).y, (*set_itr).z);

  Vec3D m_AB = (A + B).num_multi(0.5);
  Vec3D m_BC = (B + C).num_multi(0.5);
  Vec3D m_CA = (C + A).num_multi(0.5);
  cout<<C<<A<<C + A<<endl;
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
  Face m_AB_m_CA_m_BC(m_AB, m_CA, m_BC);
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
  // TODO starts here, then we go on <===========
  right.segments.insert(m_AB_B);
  right.segments.insert(m_BC_B);
  right.segments.insert(m_AB_m_CA);

  right.faces.insert(m_AB_m_BC_B);
  // ----------------------------

  Shape middle;
  middle.potins.insert(m_AB);
  middle.potins.insert(m_BC);
  middle.potins.insert(m_CA);

  midlle.segments.insert(m_AB_B);
  middle.segments.insert(m_BC_B);
  middle.segments.insert(m_AB_m_CA);

  middle.faces.insert(m_AB_m_BC_B);
  //----

  */
  int x;
  cin>>x;
  return start;
}

Shape Util::merge_shapes(Shape s, Shape v)
{
  Shape m;
  /*m.points(s.points);
  m.segments(s.segments);
  m.faces(s.faces);
  for(set<Vec3D>::iterator itr=v.points.begin(); itr != v.points.end(); itr++)
    m.points.insert((*itr));

  for(set<Segment>::iterator itr=v.segments.begin(); itr != v.segments.end(); itr++)
    m.segments.insert((*itr));

  for(set<Face>::iterator itr=v.faces.begin(); itr !=v.faces.end(); itr++)
    m.faces.insert((*itr));
  */
  return m;
}
