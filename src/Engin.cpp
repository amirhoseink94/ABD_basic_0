// openGL libraries
#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <GL/glut.h>

// my libraries
#include "../include/Engin.h"

// c++ libraries
#include <vector>
#include <cmath>
#include <algorithm>


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
    Dl_t = 1./4.;
    energy_tresh_hold = 0.05;
    distance_tresh_hold = 0.1;
  }

  vec Engin::calculate_next_Q()
  {
    // update q_mad for all of the obkects
    for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
    {
      (*itr)->update_q_mad();
    }

    vec Q_p = Q;


    vec Q_current_vec = Q;

    vec dQ_col(12*N, fill::zeros);
    //vec Q_current[N];
    //int i = 0;
    /*for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
    {
        Q_current[i] = (*itr)->q;
        i++;
    }*/

    double energy_p = 0; // calculate_energy(Q_current_vec); // previous energy
    double energy_sp = 0;
    double energy_collision = 0; //calculate_energy_collision(Q_current_vec, dQ_col); //calculating the energy of collision

    //cout<<dQ_col<<"--\n"<<energy_collision<<endl;
    energy_p += energy_collision;
    double energy = 0;
    // -----

    vec step;

    double gamma = 1;
    int counter = 0;
    int s_counter = 0;
    vec Q_temp;
    bool check_first = true;
    double old_step;
    do
    {
      energy_p = calculate_energy(Q_current_vec); // previous energy

      energy_collision = calculate_energy_collision(Q_current_vec, dQ_col); //calculating the energy of collision


      energy_p += energy_collision;

      if(check_first)
      {
        energy_sp = energy_p;
      }
      if(abs(energy_p-energy_sp)<1e-9)
        s_counter++;
      else
      {
        energy_sp = energy_p;
      }
      if(s_counter>2)
        break;
      cout<<"start the loop "<<s_counter<<endl;
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

      dQ += dQ_col; // derivatio of collision of the current Q



      // calculate the global hessian
      i = 0;
      j = 0;
      for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
      {
        H.submat(i, i, i+11, i+11) = (*itr)->E_der_der( Q_current_vec.subvec(i, i+11) );
        i = i+12;
        j++;
      }
      // Hessian at this stage will remain the same, because we chose
      // the collision potential vert carefully
      H = PSPD(H);


      step = -1 * inv(H) * dQ;
      if(check_first)
      {
        check_first = false;
        old_step = (1.0 / Dl_t) * norm(step, "inf");
      }

      cout<<"step is: "<<endl<<step<<endl;
      gamma = 1;

      Q_temp = Q_p + gamma * step;
      float toi = 10;
      cout<<"starting the search ================================"<<endl;
      toi = contact_detection_topological(Q_p, Q_temp);
      cout<<"contact detection is done: "<<toi<<endl;

      if(toi != -1)
      {
        gamma = toi;
        //Q_temp = Q_p + gamma * step;
      }

      int repeat_inner_loop=0;
      do //finding the best local place that we can go with this derivite
      {
        loop++;
        // we should go for the collision, chck for wich gamma ther is no collision
        Q_temp = Q_p + gamma * step;

        gamma = gamma * 0.25;
        energy = calculate_energy(Q_temp);
        double energy_collision_t = calculate_energy_collision(Q_temp, dQ_col); //calculating the energy of collision
        energy = energy + energy_collision_t;
        //step = step - dQ_col;
        // with respect to the new position Q_temp, we have to calculate the new energy_p
        // the overal energy(movements and orthoganility) is calculated,
        // we just need to calculate the energy of collisio based on Q_temp
        if(abs(energy - energy_p)<1e-9)
        {
          Q_current_vec = Q_temp;
          cout<<"==="<<endl<<energy<<"|"<<energy_p<<endl;
          break;
        }
        cout<<"==="<<endl<<energy<<"|"<<energy_p<<endl;
        if(abs(energy - energy_p)<1e-9)
          repeat_inner_loop++;
        if(repeat_inner_loop>=10)
          break;
      } while(energy > energy_p && loop<1024 );

      energy_p = energy;

      Q_p = Q_current_vec = Q_temp;

      double w_step = (1.0 / Dl_t) * norm(step, "inf");
      if(abs(w_step-old_step)<1e-9)
        counter++;
      else
      {
        old_step = w_step;
      }
      if(counter>2)
        break;
    } while( (1.0 / Dl_t) * norm(step, "inf") > energy_tresh_hold );
    cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<the search is over: "<<counter<<" >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
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
    // we update the each q here
    int i=0;
    for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
    {
      (*itr)->q = Q.subvec(i, i+11);
      (*itr)->apply_tranformation();
      i+=12;
    }
  }

  void Engin::apply_tranformation_temp(const vec& Q_temp)
  {
    // we update the each q here
    int i=0;
    for(auto itr = dynamic_objects.begin(); itr!=dynamic_objects.end(); itr++)
    {
      (*itr)->q = Q_temp.subvec(i, i+11);
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

    float x = min( {v0(0), f_00(0), f_01(0), f_02(0),
                    v1(0), f_10(0), f_11(0), f_12(0)});
    float y = min( {v0(1), f_00(1), f_01(1), f_02(1),
                    v1(1), f_10(1), f_11(1), f_12(1)});
    float z = min( {v0(2), f_00(2), f_01(2), f_02(2),
                    v1(2), f_10(2), f_11(2), f_12(2)});

    Eigen::Vector3<float> A(x, y, z);

    float xx = max( {v0(0), f_00(0), f_01(0), f_02(0),
                    v1(0), f_10(0), f_11(0), f_12(0)});
    float yy = max( {v0(1), f_00(1), f_01(1), f_02(1),
                    v1(1), f_10(1), f_11(1), f_12(1)});
    float zz = max( {v0(2), f_00(2), f_01(2), f_02(2),
                    v1(2), f_10(2), f_11(2), f_12(2)});
    Eigen::Vector3<float> B(xx, yy, zz);

    bool res;
    Eigen::Array3<float> err(-1, -1, -1);
    Eigen::Array3<float> err_vf = ticcd::get_numerical_error({A, B}, true, true);
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
  		err_vf, ms, toi, tolerance, t_max, max_itr, output_tolerance,
  		ticcd::DEFAULT_NO_ZERO_TOI,
  		ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);

    t = toi;
    return res;
  }

  bool Engin::contact_edge_edge(const vec& v_00, const vec& v_01, const vec& u_00, const vec& u_01,
                                  const vec& v_10, const vec& v_11, const vec& u_10, const vec& u_11, float &t)
  {
    Eigen::Vector3<float> A0(v_00(0), v_00(1), v_00(2));
  	Eigen::Vector3<float> B0(v_01(0), v_01(1), v_01(2));
  	Eigen::Vector3<float> C0(u_00(0), u_00(1), u_00(2));
  	Eigen::Vector3<float> D0(u_01(0), u_01(1), u_01(2));

    Eigen::Vector3<float> A1(v_10(0), v_10(1), v_10(2));
  	Eigen::Vector3<float> B1(v_10(0), v_10(1), v_10(2));
  	Eigen::Vector3<float> C1(u_10(0), u_10(1), u_10(2));
  	Eigen::Vector3<float> D1(u_11(0), u_11(1), u_11(2));

    float x = min( {v_00(0), v_01(0), u_00(0), u_01(0),
                    v_10(0), v_10(0), u_10(0), u_11(0)});

    float y = min( {v_00(1), v_01(1), u_00(1), u_01(1),
                    v_10(1), v_10(1), u_10(1), u_11(1)});

    float z = min( {v_00(2), v_01(2), u_00(2), u_01(2),
                    v_10(2), v_10(2), u_10(2), u_11(2)});

    Eigen::Vector3<float> A(x, y, z);

    float xx = max( {v_00(0), v_01(0), u_00(0), u_01(0),
                    v_10(0), v_10(0), u_10(0), u_11(0)});

    float yy = max( {v_00(1), v_01(1), u_00(1), u_01(1),
                    v_10(1), v_10(1), u_10(1), u_11(1)});

    float zz = max( {v_00(2), v_01(2), u_00(2), u_01(2),
                    v_10(2), v_10(2), u_10(2), u_11(2)});

    Eigen::Vector3<float> B(xx, yy, zz);

  	bool res;
    Eigen::Array3<float> err(-1, -1, -1);
    Eigen::Array3<float> err_ee = ticcd::get_numerical_error({A, B}, false, true);
    ticcd::Scalar ms = 1e-8;
    ticcd::Scalar toi;
    const ticcd::Scalar tolerance = 1e-6;
    const ticcd::Scalar t_max = 1;
    const int max_itr = 1e6;
    ticcd::Scalar output_tolerance;
    res = ticcd::edgeEdgeCCD(
  		A0, B0, C0, D0,
      A1, B1, C1, D1,
  		err_ee, ms, toi, tolerance, t_max, max_itr, output_tolerance,
  		ticcd::DEFAULT_NO_ZERO_TOI,
  		ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);

    t = toi;
    return res;
  }

  double Engin::distance_PF(vec& P, const Face& F)
  {
    vec AP = P - F.x->pos;
    vec BP = P - F.y->pos;
    vec CP = P - F.z->pos;


    vec AB = F.y->pos - F.x->pos;
    vec BC = F.z->pos - F.y->pos;
    vec CA = F.x->pos - F.z->pos;

    // check for on plan
    double v_distance = abs(dot(F.n, P - F.x->pos));
    //dQ_P = derivitives_VF_dV(F.n, P.J, v_distance, ddQ_P );
    double dist = 0;

    //mat dH = F.y->J - F.x->J;
    //mat dJ = F.z->J - F.x->J;
    //mat dK = F.z->J - F.y->J;

    //dQ_F = derivitives_VF_dF(dH, dJ, F.x->J, Q_F, P.pos, F.x->pos, v_distance, ddQ_F);

    int sgn_AB = (dot(F.n_AB, P)-dot(F.n_AB, F.x->pos) >= 0) ? 1 : 0;
    int sgn_BC = (dot(F.n_BC, P)-dot(F.n_BC, F.y->pos) >= 0) ? 1 : 0;
    int sgn_CA = (dot(F.n_CA, P)-dot(F.n_CA, F.z->pos) >= 0) ? 1 : 0;

    if(sgn_AB == 1 && sgn_BC == 1 && sgn_CA == 1)
    {
      // in this case, we do not need to do anything extra, the point is in the
      // easy position to calculate everything. S, yay!

      dist = v_distance;

    }
    else if(sgn_AB == 0 && sgn_BC == 1 && sgn_CA == 1)
    {
      // in this case, the point P is outside of the triangel, i.e., out of AB,
      // so we have to look into an extra length, and then add it fo dQ_F, and ....
      int n_A = (norm_dot(AP, AB)>0)?1:0;
      int n_B = (norm_dot(BP, -AB)>0)?1:0;


      if(n_A == 1 && n_B == 1)
      {

        double h_distance = -(dot(F.n_AB, P) - dot(F.n_AB, F.x->pos));
        dist =  sqrt(v_distance*v_distance + h_distance*h_distance);
      }
      else if(n_A == 0 && n_B == 1)
      {
        dist = norm(AP);
      }
      else if(n_A == 1 && n_B == 0)
      {
        dist = norm(BP);
      }
    }
    else if(sgn_AB == 1 && sgn_BC == 0 && sgn_CA == 1)
    {
      int n_B = (norm_dot(BP, BC)>0)?1:0;
      int n_C = (norm_dot(CP, -BC)>0)?1:0;
      if(n_B == 1 && n_C == 1)
      {
        double h_distance = -(dot(F.n_BC, P) - dot(F.n_BC, F.y->pos));
        dist = sqrt(v_distance*v_distance + h_distance*h_distance);
      }
      else if(n_B == 0 && n_C == 1)
      {
        dist = norm(BP);
      }
      else if(n_B == 1 && n_C == 0)
      {
        dist = norm(CP);
      }
    }
    else if(sgn_AB == 1 && sgn_BC == 1 && sgn_CA == 0)
    {
      int n_C = (norm_dot(CP, CA)>0)?1:0;
      int n_A = (norm_dot(AP, -CA)>0)?1:0;
      if(n_C == 1 && n_A == 1)
      {
        double h_distance = -(dot(F.n_CA, P) - dot(F.n_CA, F.z->pos));
        dist = sqrt(v_distance*v_distance + h_distance*h_distance);
      }
      else if(n_C == 0 && n_A == 1)
      {
        dist = norm(CP);
      }
      else if(n_C == 1 && n_A == 0)
      {
        dist = norm(AP);
      }
    }
    else if(sgn_AB == 1 && sgn_BC == 0 && sgn_CA == 0)
    {
      dist = norm(CP);
    }
    else if(sgn_AB == 0 && sgn_BC == 1 && sgn_CA == 0)
    {
      dist = norm(AP);
    }
    else if(sgn_AB == 0 && sgn_BC == 0 && sgn_CA == 1)
    {
      dist = norm(BP);
    }
    return dist;
  }


  double Engin::distance_VF(Particle& P, const Face& F, vec& dQ_P, vec& dQ_F)
  {
    vec z(12, fill::zeros);
    dQ_P = z;
    dQ_F = z;
    vec AP = P.pos - F.x->pos;
    vec BP = P.pos - F.y->pos;
    vec CP = P.pos - F.z->pos;


    vec AB = F.y->pos - F.x->pos;
    vec BC = F.z->pos - F.y->pos;
    vec CA = F.x->pos - F.z->pos;

    // check for on plan
    double v_distance = abs(dot(F.n, P.pos - F.x->pos));
    //dQ_P = derivitives_VF_dV(F.n, P.J, v_distance, ddQ_P );
    double dist = 0;

    //mat dH = F.y->J - F.x->J;
    //mat dJ = F.z->J - F.x->J;
    //mat dK = F.z->J - F.y->J;

    //dQ_F = derivitives_VF_dF(dH, dJ, F.x->J, Q_F, P.pos, F.x->pos, v_distance, ddQ_F);

    int sgn_AB = (dot(F.n_AB, P.pos)-dot(F.n_AB, F.x->pos) >= 0) ? 1 : 0;
    int sgn_BC = (dot(F.n_BC, P.pos)-dot(F.n_BC, F.y->pos) >= 0) ? 1 : 0;
    int sgn_CA = (dot(F.n_CA, P.pos)-dot(F.n_CA, F.z->pos) >= 0) ? 1 : 0;

    if(sgn_AB == 1 && sgn_BC == 1 && sgn_CA == 1)
    {
      // in this case, we do not need to do anything extra, the point is in the
      // easy position to calculate everything. S, yay!

      dist = v_distance;

    }
    else if(sgn_AB == 0 && sgn_BC == 1 && sgn_CA == 1)
    {
      // in this case, the point P is outside of the triangel, i.e., out of AB,
      // so we have to look into an extra length, and then add it fo dQ_F, and ....
      int n_A = (norm_dot(AP, AB)>0)?1:0;
      int n_B = (norm_dot(BP, -AB)>0)?1:0;


      if(n_A == 1 && n_B == 1)
      {

        double h_distance = -(dot(F.n_AB, P.pos) - dot(F.n_AB, F.x->pos));
        dist =  sqrt(v_distance*v_distance + h_distance*h_distance);
      }
      else if(n_A == 0 && n_B == 1)
      {
        dist = norm(AP);
      }
      else if(n_A == 1 && n_B == 0)
      {
        dist = norm(BP);
      }
    }
    else if(sgn_AB == 1 && sgn_BC == 0 && sgn_CA == 1)
    {
      int n_B = (norm_dot(BP, BC)>0)?1:0;
      int n_C = (norm_dot(CP, -BC)>0)?1:0;
      if(n_B == 1 && n_C == 1)
      {
        double h_distance = -(dot(F.n_BC, P.pos) - dot(F.n_BC, F.y->pos));
        dist = sqrt(v_distance*v_distance + h_distance*h_distance);
      }
      else if(n_B == 0 && n_C == 1)
      {
        dist = norm(BP);
      }
      else if(n_B == 1 && n_C == 0)
      {
        dist = norm(CP);
      }
    }
    else if(sgn_AB == 1 && sgn_BC == 1 && sgn_CA == 0)
    {
      int n_C = (norm_dot(CP, CA)>0)?1:0;
      int n_A = (norm_dot(AP, -CA)>0)?1:0;
      if(n_C == 1 && n_A == 1)
      {
        double h_distance = -(dot(F.n_CA, P.pos) - dot(F.n_CA, F.z->pos));
        dist = sqrt(v_distance*v_distance + h_distance*h_distance);
      }
      else if(n_C == 0 && n_A == 1)
      {
        dist = norm(CP);
      }
      else if(n_C == 1 && n_A == 0)
      {
        dist = norm(AP);
      }
    }
    else if(sgn_AB == 1 && sgn_BC == 0 && sgn_CA == 0)
    {
      dist = norm(CP);
    }
    else if(sgn_AB == 0 && sgn_BC == 1 && sgn_CA == 0)
    {
      dist = norm(AP);
    }
    else if(sgn_AB == 0 && sgn_BC == 0 && sgn_CA == 1)
    {
      dist = norm(BP);
    }

    if(dist < distance_tresh_hold)
    {

      double distance_potential = -pow(dist - distance_tresh_hold, 2) * log(dist / distance_tresh_hold);
      double distance_potential_2 = -log(dist / distance_tresh_hold);
      vec der = -log(norm(AP)) * AP + -log(norm(BP)) * BP + -log(norm(CP)) * CP;
      der = 1000 * der / norm(der);
      der = distance_potential_2 * der;
      dQ_P = P.J_T * der;
      dQ_F = -(F.x->J_T + F.y->J_T + F.z->J_T) * der;
      P.fr_value += distance_potential_2;
      //P.F_list[1] =
      //cout<<AP<<"-\n "<<BP<<"-\n "<<CP<<"-\n "<<endl;
      //cout<<"the critical vector is:"<<endl<<der<<endl;
      //cout<<"we are doomed VF "<<dist<<" "<<der<<" "<<distance_potential_2<<endl;
      return distance_potential_2;
    }
    else
    {
      P.F_list[1] = vec(3, fill::zeros);
    }
    return 0;
  }

  double Engin::distance_EE(const Segment& AB_seg, const Segment& CD_seg, vec& dQ_AB, vec& dQ_CD)
  {
    // A : AB_seg.x->pos
    // B : AB_seg.y->pos
    // C : CD_seg.x->pos
    // D : CD_seg.y->pos

    vec z(12, fill::zeros);
    dQ_AB = z;
    dQ_CD = z;

    vec AB = AB_seg.y->pos - AB_seg.x->pos;
    vec CD = CD_seg.y->pos - CD_seg.x->pos;



    vec AC = CD_seg.x->pos - AB_seg.x->pos;
    vec AD = CD_seg.y->pos - AB_seg.x->pos;
    vec BC = CD_seg.y->pos - AB_seg.x->pos;
    vec BD = CD_seg.y->pos - AB_seg.y->pos;



    mat M(2,2, fill::zeros);
    M(0,0) = norm(AB)*norm(AB);
    M(1,1) = norm(CD)*norm(CD);

    M(0,1) = dot(AB, -CD);
    M(1,0) = M(0,1);

    double dist;

    if( M(0,0)*M(1,1) - M(0,1)*M(1,0) == 0)
    {
      if(dot(AB, -AC)!=0)
      {

        double k = norm(CD)/norm(AB);
        double h = dot(AB, AC) / pow(norm(AB), 2);

        double t_1 = 0;
        double s_1 = - h / k;

        double t_2 = 1;
        double s_2 = - (h-1)/k;

        double t_3 = h;
        double s_3 = 0;

        double t_4 = h+k;
        double s_4 = 1;
        if( 0 <= s_1 && s_1 <=1)
        {
          vec P = AB_seg.x->pos + t_1 * AB;
          vec Q = CD_seg.x->pos + s_1 * CD;
          dist = norm (P - Q);

        }
        else if( 0 <= s_2 && s_2 <=1)
        {
          vec P = AB_seg.x->pos + t_2 * AB;
          vec Q = CD_seg.x->pos + s_2 * CD;
          dist = norm (P - Q);

        }
        else if( 0 <= t_3 && t_3 <=1)
        {
          vec P = AB_seg.x->pos + t_3 * AB;
          vec Q = CD_seg.x->pos + s_3 * CD;
          dist = norm (P - Q);

        }
        else if( 0 <= t_4 && t_4 <=1)
        {
          vec P = AB_seg.x->pos + t_4 * AB;
          vec Q = CD_seg.x->pos + s_4 * CD;
          dist = norm (P - Q);

        }
        else
        {
          dist = std::min( std::min( norm(AC), norm(AD) ),
          std::min( norm(BC), norm(BD) ) );

        }
      }
      else
      {
        if(norm(AB) >= norm(CD))
        {
          double t = AC(0)/AB(0);
          double s = AD(0)/AB(0);
          if ( ((0 <= t) && (t <= 1)) || ( (0 <= s) && (s <= 1)) )
          {
            dist = 0;

          }
          else
          {
            dist = std::min( std::min( norm(AC), norm(AD) ),
            std::min( norm(BC), norm(BD) ) );

          }
        }
        else
        {
          double t = -AC(0)/CD(0);
          double s = -BC(0)/CD(0);
          if ( ((0 <= t) && (t <= 1)) || ( (0 <= s) && (s <= 1)) )
          {
            dist = 0;

          }
          else
          {
            dist = std::min( std::min( norm(AC), norm(AD) ),
            std::min( norm(BC), norm(BD) ) );

          }
        }
      }
    }
    else
    {

      vec N(2, fill::zeros);
      N(0) = dot(AB, AC);
      N(1) = dot(-CD, AC);
      //cout<<M<<"\n---\n"<<det(M)<<endl;
      vec ts = pinv(M) * N;
      //cout<<"passed"<<endl;

      int sgn_AB = ( (0 < ts(0)) && (ts(0) < 1) ) ? 1 : 0;
      int sgn_CD = ( (0 < ts(1)) && (ts(1) < 1) ) ? 1 : 0;

      // case 1
      if( sgn_AB == 1 && sgn_CD == 1)
      {
        vec P = AB_seg.x->pos + ts(0) * AB;
        vec Q = CD_seg.x->pos + ts(1) * CD;

        dist = norm(P-Q);

      }
      else if(sgn_AB == 0 && sgn_CD == 1)
      {
        if( ts(0) <= 0 )
        {
          vec P = AB_seg.x->pos;
          vec Q = CD_seg.x->pos + ts(1) * CD;



          dist = norm(P-Q);

        }
        else if( 1<= ts(0) )
        {
          vec P = AB_seg.y->pos;
          vec Q = CD_seg.x->pos + ts(1) * CD;



          dist = norm(P-Q);

        }
      }
      else if(sgn_AB == 1 && sgn_CD == 0)
      {
        if( ts(1) <= 0 )
        {
          vec P = AB_seg.x->pos + ts(0) * AB;
          vec Q = CD_seg.x->pos;



          dist = norm(P-Q);

        }
        else if( 1 <= ts(1) )
        {
          vec P = AB_seg.y->pos + ts(0) * AB;
          vec Q = CD_seg.y->pos;

          dist = norm(P-Q);

        }
      }
      else if(sgn_AB == 0 && sgn_CD == 0)
      {
        vec P,Q;
        if( ts(0) <= 0 )
        {
          P = AB_seg.x->pos;
        }
        else if( 1<= ts(0) )
        {
          P = AB_seg.y->pos;
        }
        if( ts(1) <= 0 )
        {
          Q = CD_seg.x->pos;
        }
        else if( 1 <= ts(1) )
        {
          Q = CD_seg.y->pos;
        }

        dist = norm(P-Q);
      }
    }
    if(dist < distance_tresh_hold)
    {
      double distance_potential = -pow(dist - distance_tresh_hold, 2) * log(dist / distance_tresh_hold);
      double distance_potential_2 = -log(dist / distance_tresh_hold);
      vec der = -log(norm(AC)) * AC + -log(norm(AD)) * AD + -log(norm(BC)) * BC + -log(norm(BD)) * BD;
      der = 1000 * der / norm(der);
      //cout<<"we are doomed EE "<<der<<endl;
      der = distance_potential_2 * der;
      dQ_CD = (CD_seg.x->J_T + CD_seg.y->J_T ) * der;
      dQ_AB = -(AB_seg.x->J_T + AB_seg.y->J_T) * der;
      return distance_potential_2;
    }
    return 0;
  }

  double Engin::calculate_energy_collision(const vec& Q_temp, vec& dQ_col)
  {
    double energy = 0;
    dQ_col.zeros();
    apply_tranformation_temp(Q_temp);

    for(unsigned int i = 0; i<dynamic_objects.size(); i++)
    {
      // first look at all of the collisions with other dynamic_objects
      for(unsigned int j=i+1; j<dynamic_objects.size(); j++)
      {
        // each Particle of i with each face of j
        for(unsigned int i_points = 0; i_points < dynamic_objects[i]->points.size(); i_points++)
        {
          dynamic_objects[i]->points[i_points]->fr_value = 0;
          for(unsigned int j_faces = 0; j_faces < dynamic_objects[i]->faces.size(); j_faces++)
          {
            vec dQ_P, dQ_F;

            double eng = distance_VF(*dynamic_objects[i]->points[i_points], dynamic_objects[j]->faces[j_faces],
                                      dQ_P, dQ_F);
            energy = energy + eng;
            dQ_col.subvec(i, i+11) = dQ_col.subvec(i, i+11) + dQ_P;
            dQ_col.subvec(j, j+11) = dQ_col.subvec(j, j+11) + dQ_F;

          }
        } //end for on i_point

        // each segment of i with each each segment of j
        for(unsigned int i_segments = 0; i_segments < dynamic_objects[i]->segments.size(); i_segments++)
        {
          for(unsigned int j_segments = 0; j_segments < dynamic_objects[i]->segments.size(); j_segments++)
          {
            vec dQ_AB, dQ_CD;

            double eng = distance_EE(dynamic_objects[i]->segments[i_segments], dynamic_objects[j]->segments[j_segments],
                                      dQ_AB, dQ_CD);

            energy = energy + eng;
            dQ_col.subvec(i, i+11) = dQ_col.subvec(i, i+11) + dQ_AB;
            dQ_col.subvec(j, j+11) = dQ_col.subvec(j, j+11) + dQ_CD;

          }
        } // end for on i_segments
      } //end for on j for other dynamic objects

      // second, look at all of the collision with static_objects
      for(unsigned int j=0; j<static_objects.size(); j++)
      {
        // each Particle of i with each face of j
        for(unsigned int i_points = 0; i_points < dynamic_objects[i]->points.size(); i_points++)
        {
          dynamic_objects[i]->points[i_points]->fr_value = 0;
          for(unsigned int j_faces = 0; j_faces < static_objects[i]->faces.size(); j_faces++)
          {
            vec dQ_P, dQ_F;

            double eng = distance_VF(*dynamic_objects[i]->points[i_points], static_objects[j]->faces[j_faces],
                                      dQ_P, dQ_F);

            energy = energy + eng;
            dQ_col.subvec(i, i+11) = dQ_col.subvec(i, i+11) + dQ_P;
          }
        }
        //cout<<"VF is done"<<endl;
        // each segment of i with each each segment of j
        for(unsigned int i_segments = 0; i_segments < dynamic_objects[i]->segments.size(); i_segments++)
        {
          for(unsigned int j_segments = 0; j_segments < static_objects[i]->segments.size(); j_segments++)
          {
            vec dQ_AB, dQ_CD;

            double eng = distance_EE(dynamic_objects[i]->segments[i_segments], static_objects[j]->segments[j_segments],
                                      dQ_AB, dQ_CD);

            energy = energy + eng;
            dQ_col.subvec(i, i+11) = dQ_col.subvec(i, i+11) + dQ_AB;
          }
        }
      } // end for j on static objetcs
    } // end of for for all dynamic objects

    return energy;
  }
  float Engin::contact_detection_topological2(vec& Q1, vec& Q2)
  {
    vec A = dynamic_objects[0]->center->J * Q1;
    vec B = dynamic_objects[0]->center->J * Q2;
    vec P = B;
    double tresh = 0.5;

    bool not_enough = true;
    unsigned int i = 0;
    unsigned int j = 0;
    double gamma = 1;
    double gamma_0 = 0;
    double gamma_1 = 1;

    vec nc(3, fill::zeros);
    nc(1) = 75;

    float current_dist = norm(P - nc);
    if(current_dist < dynamic_objects[0]->radius)
    {
      not_enough = false;
      gamma_1 = 1;
    }
    if(not_enough)
    {
      return -1;
    }

    vec AB = B - A;
    not_enough = true;
    gamma = (gamma_0 + gamma_1)/2;
    for(int check = 1; check <10; check++)
    {
      P = A + gamma * AB;
      i = 0;
      j = 0;

      double current_dist = norm(P - nc) + tresh;
      if(current_dist < dynamic_objects[0]->radius)
      {
        not_enough = false;
        gamma_1 = gamma_0 +gamma;
      }
      if(not_enough)
      {
        gamma_0 = gamma_1 - gamma;
      }
      gamma = (gamma_0 +gamma_1)/2;
    }
    return gamma;
  }

  float Engin::contact_detection_topological(vec& Q1, vec& Q2)
  {
    vec A = dynamic_objects[0]->center->J * Q1;
    vec B = dynamic_objects[0]->center->J * Q2;
    vec P = B;
    double tresh = 0.5;

    bool not_enough = true;
    unsigned int i = 0;
    unsigned int j = 0;
    double gamma = 1;
    double gamma_0 = 0;
    double gamma_1 = 1;

    for(unsigned int j_faces = 0; j_faces < static_objects[j]->faces.size() && not_enough; j_faces++)
    {

      double current_dist = distance_PF(P, static_objects[j]->faces[j_faces]) + tresh;
      if(current_dist < dynamic_objects[0]->radius)
      {
        not_enough = false;
        gamma_1 = 1;
      }
    }
    if(not_enough)
    {
      return -1;
    }

    vec AB = B - A;
    not_enough = true;
    gamma = (gamma_0 + gamma_1)/2;
    for(int check = 1; check <10; check++)
    {
      P = A + gamma * AB;
      i = 0;
      j = 0;
      for(unsigned int j_faces = 0; j_faces < static_objects[j]->faces.size() && not_enough; j_faces++)
      {

        double current_dist = distance_PF(P, static_objects[j]->faces[j_faces]) + tresh;
        if(current_dist < dynamic_objects[0]->radius)
        {
          not_enough = false;
          gamma_1 = gamma_0 +gamma;
        }
      }
      if(not_enough)
      {
        gamma_0 = gamma_1 - gamma;
      }
      gamma = (gamma_0 +gamma_1)/2;
    }
    return gamma;
  }

  float Engin::contact_detection_exact(vec& Q1, vec& Q2)
  {
    float toi = 10;
    //cout<<"here"<<endl;
    for(unsigned int i = 0; i<dynamic_objects.size(); i++)
    {
      // first look at all of the collisions with other dynamic_objects
      for(unsigned int j=i+1; j<dynamic_objects.size(); j++)
      {
        // each Particle of i with each face of j
        for(unsigned int i_points = 0; i_points < dynamic_objects[i]->points.size(); i_points++)
        {
          for(unsigned int j_faces = 0; j_faces < dynamic_objects[j]->faces.size(); j_faces++)
          {
            vec v1 = dynamic_objects[i]->points[i_points]->J * Q1.subvec(i, i+11);

            vec f11 = dynamic_objects[j]->faces[j_faces].x->J * Q1.subvec(j, j+11);
            vec f12 = dynamic_objects[j]->faces[j_faces].y->J * Q1.subvec(j, j+11);
            vec f13 = dynamic_objects[j]->faces[j_faces].z->J * Q1.subvec(j, j+11);

            vec v2 = dynamic_objects[i]->points[i_points]->J * Q2.subvec(i, i+11);

            vec f21 = dynamic_objects[j]->faces[j_faces].x->J * Q2.subvec(j, j+11);
            vec f22 = dynamic_objects[j]->faces[j_faces].y->J * Q2.subvec(j, j+11);
            vec f23 = dynamic_objects[j]->faces[j_faces].z->J * Q2.subvec(j, j+11);
            float temp_t;

            bool res_contact = contact_vertex_face(v1, f11, f12, f13, v2, f21, f22, f23, temp_t);
            if(res_contact && temp_t < toi)
            {
              toi = temp_t;
            }
          }
        } //end for on i_point

        // each segment of i with each each segment of j
        for(unsigned int i_segments = 0; i_segments < dynamic_objects[i]->segments.size(); i_segments++)
        {
          for(unsigned int j_segments = 0; j_segments < dynamic_objects[j]->segments.size(); j_segments++)
          {
            vec A0 = dynamic_objects[i]->segments[i_segments].x->J * Q1.subvec(i, i+11);
            vec B0 = dynamic_objects[i]->segments[i_segments].y->J * Q1.subvec(i, i+11);
            vec C0 = dynamic_objects[j]->segments[j_segments].x->J * Q1.subvec(j, j+11);
            vec D0 = dynamic_objects[j]->segments[j_segments].y->J * Q1.subvec(j, j+11);

            vec A1 = dynamic_objects[i]->segments[i_segments].x->J * Q2.subvec(i, i+11);
            vec B1 = dynamic_objects[i]->segments[i_segments].y->J * Q2.subvec(i, i+11);
            vec C1 = dynamic_objects[j]->segments[j_segments].x->J * Q2.subvec(j, j+11);
            vec D1 = dynamic_objects[j]->segments[j_segments].y->J * Q2.subvec(j, j+11);
            float temp_t;
            bool res_contact = contact_vertex_face(A0, B0, C0, D0, A1, B1, C1, D1, temp_t);
            if(res_contact && temp_t < toi)
            {
              toi = temp_t;
            }

          }
        } // end for on i_segments
      } //end for on j for other dynamic objects

      // second, look at all of the collision with static_objects
      for(unsigned int j=0; j<static_objects.size(); j++)
      {
        // each Particle of i with each face of j
        for(unsigned int i_points = 0; i_points < dynamic_objects[i]->points.size(); i_points++)
        {
          for(unsigned int j_faces = 0; j_faces < static_objects[j]->faces.size(); j_faces++)
          {

            vec v1 = dynamic_objects[i]->points[i_points]->J * Q1.subvec(i, i+11);

            vec f11 = static_objects[j]->faces[j_faces].x->pos;
            vec f12 = static_objects[j]->faces[j_faces].y->pos;
            vec f13 = static_objects[j]->faces[j_faces].z->pos;

            vec v2 = dynamic_objects[i]->points[i_points]->J * Q2.subvec(i, i+11);

            float temp_t;
            bool res_contact = contact_vertex_face(v1, f11, f12, f13, v2, f11, f12, f13, temp_t);
            if(res_contact && temp_t < toi)
            {
              toi = temp_t;
            }
          }
        }

        // each segment of i with each each segment of j
        for(unsigned int i_segments = 0; i_segments < dynamic_objects[i]->segments.size(); i_segments++)
        {
          for(unsigned int j_segments = 0; j_segments < static_objects[j]->segments.size(); j_segments++)
          {
            vec A0 = dynamic_objects[i]->segments[i_segments].x->J * Q1.subvec(i, i+11);
            vec B0 = dynamic_objects[i]->segments[i_segments].y->J * Q1.subvec(i, i+11);
            vec C0 = static_objects[j]->segments[j_segments].x->pos;
            vec D0 = static_objects[j]->segments[j_segments].y->pos;

            vec A1 = dynamic_objects[i]->segments[i_segments].x->J * Q2.subvec(i, i+11);
            vec B1 = dynamic_objects[i]->segments[i_segments].y->J * Q2.subvec(i, i+11);

            float temp_t;
            //cout<<i_segments<< " " <<j_segments<<endl;
            bool res_contact = contact_vertex_face(A0, B0, C0, D0, A1, B1, C0, D0, temp_t);
            if(res_contact && temp_t < toi)
            {
              toi = temp_t;
            }
          }
        }
      } // end for j on static objetcs
    } // end of for for all dynamic objects

    return toi;
  }

}
