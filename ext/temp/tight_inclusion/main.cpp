#include <stdio.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <vector>
// my libraries

#include "ccd.hpp"
#include "timer.hpp"
//#include "rational/ccd.hpp"
//#include "interval_root_finder.hpp"
//namespaces
using namespace ticcd;
using namespace std;


int main( int argc, char** argv)
{
	
	Eigen::Vector3<float> v(1,1,3);

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
	//while(v2.x() > -1)
	{
		ticcd::Scalar output_tolerance;
		res = ticcd::vertexFaceCCD(
				v,
				f00, f01, f02,
				v2,
				f00, f01, f02,
				err, ms, toi, tolerance, t_max, max_itr, output_tolerance,
				ticcd::DEFAULT_NO_ZERO_TOI,
				ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);

		if(res)
		{
			cout<<"we have a contact"<<endl;
			cout<<"time of impact:"<< toi <<endl;
			int x;
			cin>>x;
		}
		else
			cout<<"we are good, keep travelling!"<<endl;
		cout<< "the movement is toward: "<< v2 <<endl;
		//v2.z() = v2.z() - 0.01;
	}
	cout<<"the program ends"<<endl;
	return 0;
}
