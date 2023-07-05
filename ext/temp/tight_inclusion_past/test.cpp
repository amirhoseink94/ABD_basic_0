#include <iostream>
//#include <armadillo>
#include <cmath>
#include <tight_inclusion/ccd.hpp>
//using namespace arma;
using namespace std;
#define PI 3.14159265
int main()
{
	//mat A(5, 5, fill::randu);
	//cout<<A<<endl<<"==="<<endl;
	cout<<"thisgs are good!"<<endl;
	Eigen::Vector3<float> v(1,1,1);

	Eigen::Vector3<float> f00(0,0,0);
	Eigen::Vector3<float> f01(10,0,0);
	Eigen::Vector3<float> f02(0,10,0);

	Eigen::Vector3<float> v2(1,1,-1);
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
		v2,
		f00, f01, f02,
		err, ms, toi, tolerance, t_max, max_itr, output_tolerance,
		ticcd::DEFAULT_NO_ZERO_TOI,
		ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);

	cout<<"res "<<res<<endl;
	cout<<"the program starts"<<endl;
	//A = 10 * 10  *A;  // generate a symmetric matrix
	
	/*B(0,0) = -10;
	vec eigval;
	mat eigvec;
	

	eig_sym(eigval, eigvec, B);
	cout<< eigval << "===" <<endl;
	mat C = eigvec * diagmat(eigval) * eigvec.t();
	cout<<C<<endl;	
	for(int i=0;i<eigval.size();i++)
	{
		if(eigval(i)<0)
			eigval(i)= 0.0;
	}*/
	
	//cout<<A<<endl<<"==="<<endl;
	
	//vec eigvalp;
	//mat eigvecp;
	

	//eig_sym(eigvalp, eigvecp, C);	
	//cout<<eigvalp<<"==="<<eigvecp<<endl;
	
	return 0;
}

