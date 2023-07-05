#include "../include/Utilities.h"

#include <armadillo>

using namespace std;
using namespace arma;

typedef Mat<double> Vec3f;

Body Util::make_mesh(Body start, int n)
{
    static Body res = Util::make_mesh_handler(start, 0, n);
    int x;
    cout<<"the moment of truth"<<endl;
    cin>>x;
    cout<<res.segments[0];
    return res;
}

Body Util::make_mesh_handler(Body& start, int itr, int n)
{
  cout<<"itr: "<<itr<<endl;
  if(itr == n)
  {
    return start;
  }
  cout<<"we are here in the handler"<<endl;
  Vec3f Z(3,1, fill::zeros);
  /*static Particle A(start.points[0].pos, Z, 1, Z);
  static Particle B(start.points[1].pos, Z, 1, Z);
  static Particle C(start.points[2].pos, Z, 1, Z);
  Vec3f p1 = ((A.pos + B.pos)*(0.5));
  Vec3f p2 = ((B.pos + C.pos)*(0.5));
  Vec3f p3 = ((C.pos + A.pos)*(0.5));
  static Particle m_AB(p1, Z, 1, Z);
  static Particle m_BC(p2, Z, 1, Z);
  static Particle m_CA(p3, Z, 1, Z);
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

  cout<<"the up body"<<endl;
  up.print();
  int xx;
  cin>>xx;
  // ----------------------------

  Body left;
  left.points.push_back(B);
  left.points.push_back(m_AB);
  left.points.push_back(m_BC);

  left.segments.push_back(m_AB_B);
  left.segments.push_back(m_BC_B);
  left.segments.push_back(m_AB_m_CA);

  left.faces.push_back(m_AB_m_BC_B);

  cout<<"the left body"<<endl;
  left.print();

  cin>>xx;
  // ----------------------------

  Body right;
  right.points.push_back(C);
  right.points.push_back(m_BC);
  right.points.push_back(m_CA);

  right.segments.push_back(m_BC_C);
  right.segments.push_back(m_CA_C);
  right.segments.push_back(m_BC_m_CA);

  right.faces.push_back(m_BC_m_CA_C);
  cout<<"the right body"<<endl;
  right.print();

  cin>>xx;
  // ----------------------------

  Body middle;
  middle.points.push_back(m_AB);
  middle.points.push_back(m_BC);
  middle.points.push_back(m_CA);

  middle.segments.push_back(m_AB_m_CA);
  middle.segments.push_back(m_BC_m_CA);
  middle.segments.push_back(m_AB_m_BC);

  middle.faces.push_back(m_AB_m_BC_m_CA);
  cout<<"the middle body"<<endl;
  middle.print();

  cin>>xx;
  //----
  Body res_up = make_mesh_handler(up, itr+1, n);
  Body res_middle = make_mesh_handler(middle, itr+1, n);
  Body res_right = make_mesh_handler(right, itr+1, n);
  Body res_left = make_mesh_handler(left, itr+1, n);
  static Body res = merge_shapes(merge_shapes(res_left, res_right), merge_shapes(res_up, res_middle));


  //res.update();
  cout<<"we made out of update!"<<endl;
  cout<<"<><><><><><><>"<<endl;
  res.print();
  cout<<res.segments.size()<<"|"<<res.faces.size()<<endl;
  cout<<res.segments[0]<<endl;
  cout<<"<><><><><><><>DONE"<<endl;
  int xxx;
  cin>>xxx;*/
  return start;
}

Body Util::merge_shapes(Body A, Body B)
{
  Body m;
  m.points.insert(m.points.begin(), A.points.begin(), A.points.end());
  m.segments.insert(m.segments.begin(), A.segments.begin(), A.segments.end());
  m.faces.insert(m.faces.begin(), A.faces.begin(), A.faces.end());
  cout<<"just this time"<<endl;
  cout<<m.segments[0]<<endl;
  m.points.insert(m.points.end(), B.points.begin(), B.points.end());
  m.segments.insert(m.segments.end(), B.segments.begin(), B.segments.end());
  m.faces.insert(m.faces.end(), B.faces.begin(), B.faces.end());

  return m;
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
