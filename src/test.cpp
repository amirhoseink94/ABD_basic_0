#include <iostream>
#include <armadillo>
#include <cmath>
#include "../include/Segment.h"
#include "../include/Particle.h"

using namespace arma;
using namespace std;
#define PI 3.14159265
int main()
{
	Mat<float> A(3,1, fill::zeros);
	Mat<float> B(3,1, fill::zeros);
	Mat<float> Z(3,1, fill::zeros);
	A[0] = 1;
	A[1] = 2;
	A[2] = 3;
	
	Mat<float> C = A+2;
	Particle a(A, Z, 10, Z);
	Particle b(B, Z, 10, Z);
	Particle c(C, Z, 10, Z);
	Segment t(&a, &b);
	Face f(&a, &b, &c);
	cout<<a<<endl<<b<<endl;
	cout<<t<<endl;
	a.pos = a.pos+10;
	cout<<"next"<<endl;
	cout<<t<<endl;
	
	cout<<f<<endl;
	
	return 0;
}

