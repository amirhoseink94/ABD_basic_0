#include "../include/Face.h"

using namespace std;

Face::Face(Particle* x, Particle* y, Particle* z)
{
  this->x = x;
  this->y = y;
  this->z = z;
  C_0 = vec(3, fill::zeros);
  c_0 = vec(3, fill::zeros);
  J_C = mat(3, 12, fill::zeros);
  J_c = mat(3, 12, fill::zeros);
  n = vec(3, fill::zeros);
  n_AB = vec(3, fill::zeros);
  n_CA = vec(3, fill::zeros);

  Cal_centre();
}

void Face::Cal_centre()
{
  n = cross(y->pos - x->pos, z->pos - x->pos);
  n = n / norm(n);
  vec M = 0.5*(x->pos + y->pos);
  vec N = 0.5*(x->pos + z->pos);

  n_AB = -cross(y->pos - x->pos, n);
  n_CA = cross(z->pos - x->pos, n);
  n_AB = n_AB / norm(n_AB);
  n_CA = n_CA / norm(n_CA);

  n_BC = -cross(z->pos - y->pos, n);
  n_BC = n_BC / norm(n_BC);



  /*mat G(2,2, fill::zeros);
  G(0,0) = n_AB[0];
  G(1,0) = n_AB[1];
  G(0,1) = -n_CA[0];
  G(1,1) = -n_CA[1];

  vec T(2, fill::zeros);
  T(0) = N[0] - M[0];
  T(1) = N[1] - M[1];

  vec ts = inv(G) * T;



  C_0[0] = ts[0] * n_AB[0] + M[0];
  C_0[1] = ts[0] * n_AB[1] + M[1];
  C_0[2] = ts[0] * n_AB[2] + M[2];
  C = C_0;
  R = norm(C_0 - x->pos);

  J_C(0,0) = 1;
	J_C(1,1) = 1;
	J_C(2,2) = 1;


	J_C(0,3) = C_0(0);
	J_C(1,4) = C_0(0);
	J_C(2,5) = C_0(0);

	J_C(0,6) = C_0(1);
	J_C(1,7) = C_0(1);
	J_C(2,8) = C_0(1);

	J_C(0,9) = C_0(2);
	J_C(1,10) = C_0(2);
	J_C(2,11) = C_0(2);


  vec BC = z->pos - y->pos;
  BC = BC/norm(BC);
  vec CB = -BC;
  vec BA = x->pos - y->pos;
  BA = BA/norm(BA);
  vec CA = x->pos - z->pos;
  CA = CA/norm(CA);

  vec t1 = (BC + BA) * 0.5;
  vec t2 = (CB + CA) * 0.5;

  G(0,0) = t2[0];
  G(1,0) = t2[1];
  G(0,1) = -t1[0];
  G(1,1) = -t1[1];

  T(0) = y->pos[0] - z->pos[0];
  T(1) = y->pos[1] - z->pos[1];

  ts = inv(G) * T;
  c_0[0] = ts[1] * t1[0] + y->pos[0];
  c_0[1] = ts[1] * t1[1] + y->pos[1];
  c_0[2] = ts[1] * t1[2] + y->pos[2];

  c = c_0;
  J_c(0,0) = 1;
	J_c(1,1) = 1;
	J_c(2,2) = 1;


	J_c(0,3) = c_0(0);
	J_c(1,4) = c_0(0);
	J_c(2,5) = c_0(0);

	J_c(0,6) = c_0(1);
	J_c(1,7) = c_0(1);
	J_c(2,8) = c_0(1);

	J_c(0,9) = c_0(2);
	J_c(1,10) = c_0(2);
	J_c(2,11) = c_0(2);

  r = norm(cross(c_0 - y->pos, BA)) / norm(BA);*/
  return;
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
