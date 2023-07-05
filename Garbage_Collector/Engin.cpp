// openGL libraries
#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <GL/glut.h>

// my libraries
#include "../include/Engin.h"

// c++ libraries
#include <vector>
#include <cmath>

// namespaces
using namespace std;
using namespace arma;

namespace ABD
{
  void Engin::add_object(Body* b)
  {
    if(b->movable)
      dynamic_objects.push_back(b);
    else
      static_objects.push_back(b);
  }

  void Engin::init()
  {
    N = dynamic_objects.size();

    Q.zeros(12 * N);
    dQ.zeros(12 * N);
    H.zeros(12 * N, 12 * N);
    int i = 0;
    for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
    {
        Q.subvec(i, i+11) = (*itr)->q;
        i = i+12;
    }
    Dl_t = 1./32.;
    energy_tresh_hold = 0.05;
  }

  vec Engin::calculate_next_Q()
  {
    // update q_mad for all of the obkects
    for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
    {
      (*itr)->update_q_mad();
    }

    mat Q_p = Q;


    vec Q_current_vec = Q;
    //vec Q_current[N];
    int i = 0;
    /*for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
    {
        Q_current[i] = (*itr)->q;
        i++;
    }*/

    double energy_p = calculate_energy(Q_current_vec); // previous energy
    double energy = 0;
    // -----

    vec step;

    double gamma = 1;
    int counter = 0;

    Mat<double> Q_temp;

    do
    {
      int loop = 0;
      //calculate te global derivation
      int i = 0;
      int j = 0;
      for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
      {
        dQ.subvec(i, i+11) = (*itr)->E_der( Q_current_vec.subvec(i, i+11) );
        i = i+12;
        j++;
      }

      // calculate the global hessian
      i = 0;
      j = 0;
      for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
      {
        H.submat(i, i, i+11, i+11) = (*itr)->E_der_der( Q_current_vec.subvec(i, i+11) );
        i = i+12;
        j++;
      }

      H = PSPD(H);

      step = -1 * inv(H) * dQ;

      cout<<"step is: "<<endl<<step<<endl;

      gamma = 1;

      do
      {
        loop++;
        Q_temp = Q_p + gamma * step;
        //q_current = q_p + gamma * step;

        gamma = gamma * 0.5;
        energy = calculate_energy(Q_temp);
        if(energy<energy_p)
        {
          Q_current_vec = Q_temp;
          cout<<"==="<<endl<<energy<<"|"<<energy_p<<endl;
          break;
        }
        cout<<"==="<<endl<<energy<<"|"<<energy_p<<endl;

      } while(energy > energy_p && loop<1024 );

      energy_p = energy;

      Q_p = Q_current_vec;

      counter++;
      if(counter>1000)
        break;

    } while((1.0 / Dl_t) * norm(step, "inf")  > energy_tresh_hold );

    return Q_current_vec;

  }


  mat Engin::PSPD(const mat& X)
  {
    vec eigval;
  	mat eigvec;

  	eig_sym(eigval, eigvec, X);
    //cout<<"the init is:"<<endl<<eigval<<endl;
    for(unsigned long int i=0;i<eigval.size();i++)
  	{
  		if(eigval(i)<0)
  			eigval(i) = 0.000001;
  	}

  	mat C = eigvec * diagmat(eigval) * eigvec.t();
    return C;
  }
  double Engin::calculate_energy(const vec& X)
  {
    double global_e = 0;
    int i = 0;
    for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
    {
      global_e += (*itr)->calculate_energy(X.subvec(i, i+11));
      i = i+12;
    }
    return global_e;
  }
  void Engin::apply_tranformation()
  {
    int i=0;
    for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
    {
      (*itr)->q = Q.subvec(i, i+11);
      (*itr)->apply_tranformation();
      i+=12;
    }
  }

  bool Engin::contact_vertex_face(const vec& v0, const vec& f_00, const vec& f_01, const vec& f_02,
                                  const vec& v1, const vec& f_10, const vec& f_11, const vec& f_12, float &t)
  {
    Eigen::Vector3<float> v(v0(0),v0(1),v0(2));

  	Eigen::Vector3<float> f00(f_00(0), f_00(1), f_00(2));
  	Eigen::Vector3<float> f01(f_01(0), f_01(1), f_01(2));
  	Eigen::Vector3<float> f02(f_02(0), f_02(1), f_02(2));

  	Eigen::Vector3<float> vp(v1(0),v1(1),v1(2));

    Eigen::Vector3<float> f10(f_10(0), f_10(1), f_10(2));
  	Eigen::Vector3<float> f11(f_11(0), f_11(1), f_11(2));
  	Eigen::Vector3<float> f12(f_12(0), f_12(1), f_12(2));

  	bool res;
    Eigen::Array3<float> err(-1, -1, -1);
    ticcd::Scalar ms = 1e-8;
    ticcd::Scalar toi;
    const ticcd::Scalar tolerance = 1e-6;
    const ticcd::Scalar t_max = 1;
    const int max_itr = 1e6;
    ticcd::Scalar output_tolerance;
    res = ticcd::vertexFaceCCD(
  		v,
  		f00, f01, f02,
  		vp,
  		f00, f01, f02,
  		err, ms, toi, tolerance, t_max, max_itr, output_tolerance,
  		ticcd::DEFAULT_NO_ZERO_TOI,
  		ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);

    t = toi;
    return;
  }


  double Engin::distance_VF(const Particle& P, const Face& F, vec& dQ_P, vec& dQ_F
                                                            , mat& ddQ_P, mat& ddQ_F
                                                            , const vec& Q_P, const vec& Q_F)
  {
    vec AP = P.pos - F.x->pos;
    vec BP = P.pos - F.y->pos;
    vec CP = P.pos - F.z->pos;

    cout<<P.pos<<endl;
    cout<<F.x->pos<<endl<<F.y->pos<<endl<<F.z->pos;

    vec AB = F.y->pos - F.x->pos;
    vec BC = F.z->pos - F.y->pos;
    vec CA = F.x->pos - F.z->pos;

    // check for on plan
    double v_distance = (dot(F.n, P.pos) - dot(F.n, F.x->pos));
    dQ_P = derivitives_VF_dV(F.n, P.J, v_distance, ddQ_P );

    mat dH = F.y->J - F.x->J;
    mat dJ = F.z->J - F.x->J;
    mat dK = F.z->J - F.y->J;

    dQ_F = derivitives_VF_dF(dH, dJ, F.x->J, Q_F, P.pos, F.x->pos, v_distance, ddQ_F);

    int sgn_AB = (dot(F.n_AB, P.pos)-dot(F.n_AB, F.x->pos) >= 0) ? 1 : 0;
    int sgn_BC = (dot(F.n_BC, P.pos)-dot(F.n_BC, F.y->pos) >= 0) ? 1 : 0;
    int sgn_CA = (dot(F.n_CA, P.pos)-dot(F.n_CA, F.z->pos) >= 0) ? 1 : 0;

    if(sgn_AB == 1 && sgn_BC == 1 && sgn_CA == 1)
    {
      // in this case, we do not need to do anything extra, the point is in the
      // easy position to calculate everything. S, yay!
      cout<<"the easy part!"<<endl;
      double d = v_distance* v_distance;
      return d;
    }
    else if(sgn_AB == 0 && sgn_BC == 1 && sgn_CA == 1)
    {
      // in this case, the point P is outside of the triangel, i.e., out of AB,
      // so we have to look into an extra length, and then add it fo dQ_F, and ....
      int n_A = (norm_dot(AP, AB)>0)?1:0;
      int n_B = (norm_dot(BP, -AB)>0)?1:0;
      cout<<"the most tricky part!"<<endl;
      if(n_A == 1 && n_B == 1)
      {
        cout<<"the most most tricky part!"<<endl;
        double h_distance = -(dot(F.n_AB, P.pos) - dot(F.n_AB, F.x->pos));
        double dist = dist_PE_P(P, *F.x, *F.y, dQ_P, ddQ_P);
        double dist2 = dist_PE_E(P, *F.x, *F.y, dQ_F, ddQ_F);
        return (v_distance*v_distance + h_distance*h_distance);
      }
      else if(n_A == 0 && n_B == 1)
      {
        double dist = dist_PP(F.x->J, F.x->pos
                            , P.J, P.pos
                            , dQ_F, ddQ_F
                            , dQ_P, ddQ_P);
        //return norm(AP)*norm(AP);
        return dist;
      }
      else if(n_A == 1 && n_B == 0)
      {
        double dist = dist_PP(F.y->J, F.y->pos
                            , P.J, P.pos
                            , dQ_F, ddQ_F
                            , dQ_P, ddQ_P);
        return norm(BP)*norm(BP);
      }
    }
    else if(sgn_AB == 1 && sgn_BC == 0 && sgn_CA == 1)
    {
      int n_B = (norm_dot(BP, BC)>0)?1:0;
      int n_C = (norm_dot(CP, -BC)>0)?1:0;
      if(n_B == 1 && n_C == 1)
      {
        double h_distance = -(dot(F.n_BC, P.pos) - dot(F.n_BC, F.y->pos));
        double dist = dist_PE_P(P, *F.y, *F.z, dQ_P, ddQ_P);
        double dist2 = dist_PE_E(P, *F.y, *F.z, dQ_F, ddQ_F);
        return (v_distance*v_distance + h_distance*h_distance);
      }
      else if(n_B == 0 && n_C == 1)
      {
        double dist = dist_PP(F.y->J, F.y->pos
                            , P.J, P.pos
                            , dQ_F, ddQ_F
                            , dQ_P, ddQ_P);
        return norm(BP)*norm(BP);
      }
      else if(n_B == 1 && n_C == 0)
      {
        double dist = dist_PP(F.z->J, F.z->pos
                            , P.J, P.pos
                            , dQ_F, ddQ_F
                            , dQ_P, ddQ_P);
        return norm(CP)*norm(CP);
      }
    }
    else if(sgn_AB == 1 && sgn_BC == 1 && sgn_CA == 0)
    {
      int n_C = (norm_dot(CP, CA)>0)?1:0;
      int n_A = (norm_dot(AP, -CA)>0)?1:0;
      if(n_C == 1 && n_A == 1)
      {
        double h_distance = -(dot(F.n_CA, P.pos) - dot(F.n_CA, F.z->pos));
        double dist = dist_PE_P(P, *F.z, *F.x, dQ_P, ddQ_P);
        double dist2 = dist_PE_E(P, *F.z, *F.x, dQ_F, ddQ_F);
        return (v_distance*v_distance + h_distance*h_distance);
      }
      else if(n_C == 0 && n_A == 1)
      {
        double dist = dist_PP(F.z->J, F.z->pos
                            , P.J, P.pos
                            , dQ_F, ddQ_F
                            , dQ_P, ddQ_P);
        return norm(CP)*norm(CP);
      }
      else if(n_C == 1 && n_A == 0)
      {
        double dist = dist_PP(F.x->J, F.x->pos
                            , P.J, P.pos
                            , dQ_F, ddQ_F
                            , dQ_P, ddQ_P);
        return norm(AP)*norm(AP);
      }
    }
    else if(sgn_AB == 1 && sgn_BC == 0 && sgn_CA == 0)
    {
      double dist = dist_PP(F.z->J, F.z->pos
                          , P.J, P.pos
                          , dQ_F, ddQ_F
                          , dQ_P, ddQ_P);
      return norm(CP)*norm(CP);
    }
    else if(sgn_AB == 0 && sgn_BC == 1 && sgn_CA == 0)
    {
      double dist = dist_PP(F.x->J, F.x->pos
                          , P.J, P.pos
                          , dQ_F, ddQ_F
                          , dQ_P, ddQ_P);
      return norm(AP)*norm(AP);
    }
    else if(sgn_AB == 0 && sgn_BC == 0 && sgn_CA == 1)
    {
      double dist = dist_PP(F.y->J, F.y->pos
                          , P.J, P.pos
                          , dQ_F, ddQ_F
                          , dQ_P, ddQ_P);
      return norm(BP)*norm(BP);
    }
    return -1;
  }
  vec Engin::derivitives_VF_dV(const vec& n, const mat& J, const double dist, mat& ddQ_P)
  {
    vec dQ_P = (2 * dist * (n.t() * J) ).t();
    ddQ_P = 2 * (n.t() * J).t() * (n.t() * J);
    return dQ_P;
  }
  vec Engin::derivitives_VF_dF(const mat& H, const mat& J, const mat& K, const vec& q, const vec& P, const vec& A, const double dist, mat& ddQ_P)
  {
    vec dQ_P = (2 * dist) * vec(12, fill::ones);
    double len = 1.0/norm(cross(H *q, J * q));
    cout<<"so far so good"<<endl;
    for(int i=0;i<12;i++)
    {
      dQ_P(i) = len * dQ_P(i) *( dot(cross(H * q , J * q), -K.col(i)) + dot(cross(H.col(i),J*q) + cross(H*q , J.col(i)) ,P-A) );
    }
    ddQ_P = dQ_P * dQ_P.t();
    for(int i=0;i<12; i++)
      for(int j=0;j<12;j++)
      {
        ddQ_P(i,j) += dot(cross(H.col(j),J*q) + cross(H*q , J.col(j)) ,-K.col(i)) +
                    dot(cross(H.col(i), J.col(j)) + cross(H.col(j),J.col(i)) , P-A) +
                    dot(cross(H.col(i),J*q) + cross(H*q , J.col(i)) ,-K.col(j));
      }
    ddQ_P = 2 * len * ddQ_P;
    return dQ_P;
  }
  double Engin::dist_PP(const mat& J_1, const vec& pos_1
                      , const mat& J_2, const vec& pos_2
                      , vec& dq_1, mat& ddq_1
                      , vec& dq_2, mat& ddq_2)
  {
    double dist = norm(pos_1 - pos_2);
    dq_1 = (dist)*(2 * (pos_1 - pos_2).t() * J_1).t();
    dq_2 = (dist)*(2 * (pos_2 - pos_1).t() * J_2).t();
    ddq_1 = J_1.t() * J_1;
    ddq_2 = J_2.t() * J_2;
    dist *= dist;
    return dist;
  }
  double Engin::dist_PE_P(const Particle& P, const Particle& A, const Particle& B
                 , vec& dq_P, mat& ddq_P)
  {
    vec b = B.pos - A.pos;
    b = b / norm(b);
    vec P_A = P.pos - A.pos;
    double dist = norm(cross(P_A, b));
    dq_P = (2 * cross( b, cross(P_A, b) ).t() * P.J).t();
    for(int i=0;i<12;i++)
    {
      vec Fp(3, fill::zeros);
      for(int j=0;j<12;j++)
      {
        Fp = 2 * cross(b, cross(P.J.col(j), b));
        ddq_P(i,j) = dot(Fp, P.J.col(i));
      }
    }
    return dist * dist;
  }
  double Engin::dist_PE_E(const Particle& P, const Particle& A, const Particle& B
                 , vec& dq_E, mat& ddq_E)
  {
    vec a = P.pos - A.pos;
    vec BA = B.pos - A.pos;
    double len = norm(BA);

    BA = (1/len)*BA;
    mat J = B.J - A.J;
    mat Fp(3, 12, fill::zeros);

    for(int i=0; i<12;i++)
    {
      Fp.col(i) = (1/len) * ( (cross( -A.J.col(i) , B.pos - A.pos ) + cross((P.pos - A.pos), J.col(i))));
    }

    double dist = norm(cross(a, BA));
    cout<<"this is important:\n"<< 2 * cross( a, cross(a, BA))<<endl;
    dq_E = (2 * (cross(a, BA)).t() * Fp).t();
    for(int i=0;i<12;i++)
    {
      for(int j=0;j<12;j++)
      {
        ddq_E(i,j) = dot((1/len) * ( (cross( -A.J.col(j) , B.pos - A.pos ) + cross((P.pos - A.pos), J.col(j)))), Fp.col(i)) +
         dot(cross(a, BA) , (1/len) *(cross(-A.J.col(i), J.col(j)) + cross(-A.J.col(j), J.col(i))));
      }
    }
    return dist * dist;
 }
}
