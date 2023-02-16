#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <cmath>
using namespace std;

template<typename T>
class Vector
{
	public:
    T x, y, z;
    Vector();
		//Vector(const Vector<T>&);
    Vector(T x, T y, T z);
    // Overload the * operator
    Vector<T> operator+ (const Vector<T>& obj);
    Vector<T> operator- (const Vector<T>& obj);
    Vector<T>& operator= (const Vector<T>& obj);
    T& operator[](int);

    float length();
    Vector<T> num_multi(float);

    template <typename U>
    friend ostream& operator<<(ostream& os, const Vector<U>& obj);
    template <typename U>
    friend bool operator== (const Vector<U>& v1, const Vector<U>& v2);
    template <typename U>
    friend bool operator!= (const Vector<U>& v1, const Vector<U>& v2);
    template <typename U>
    friend bool operator<= (const Vector<U>& v1, const Vector<U>& v2);
    template <typename U>
    friend bool operator>= (const Vector<U>& v1, const Vector<U>& v2);
    template <typename U>
    friend bool operator< (const Vector<U>& v1, const Vector<U>& v2);
    template <typename U>
    friend bool operator> (const Vector<U>& v1, const Vector<U>& v2);

};

typedef Vector<float> Vec3D;

template<typename T>
Vector<T>::Vector()
{
	this->x = 0;
	this->y = 0;
	this->z = 0;
}
/*template<typename T>
Vector<T>::Vector(const Vector<T>& obj)
{
	this->x = obj.x;
	this->y = obj.y;
	this->z = obj.z;
}*/
template<typename T>
Vector<T>::Vector(T x, T y, T z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}
template<typename T>
Vector<T> Vector<T>::operator+(const Vector<T>& obj)
{
    Vector<T> temp;
    temp.x = x + obj.x;
    temp.y = y + obj.y;
    temp.z = z + obj.z;
    return temp;
}
template<typename T>
Vector<T> Vector<T>::operator-(const Vector<T>& obj)
{
    Vector<T> temp;
    temp.x = x - obj.x;
    temp.y = y - obj.y;
    temp.z = z - obj.z;
    return temp;
}
template<typename T>
Vector<T>& Vector<T>::operator= (const Vector<T>& obj)
{
	x = obj.x;
	y = obj.y;
	z = obj.z;
  return *this;
}
template<typename T>
ostream& operator<<(ostream& os, const Vector<T>& obj)
{
    os << "(" << obj.x << ", " << obj.y << ", " << obj.z << ")";;
    return os;
}
template<typename T>
Vector<T> Vector<T>::num_multi(float a)
{
	Vector<T> temp;
    temp.x = a * x;
    temp.y = a * y;
    temp.z = a * z;

    return temp;
}
template<typename T>
float Vector<T>::length()
{
	float l = pow(x, 2) + pow(y, 2) + pow(z, 2);
	l = sqrt(l);
	return l;
}
template<typename T>
bool operator== (const Vector<T>& v1, const Vector<T>& v2)
{
	if(v1.x == v2.x && v1.y == v2.y && v1.z == v2.z)
		return true;
	else
		return false;
}
template<typename T>
bool operator!= (const Vector<T>& v1, const Vector<T>& v2)
{
	if(v1.x != v2.x || v1.y != v2.y || v1.z != v2.z)
		return true;
	else
		return false;
}
template<typename T>
bool operator<= (const Vector<T>& v1, const Vector<T>& v2)
{
	if(v1.x <= v2.x && v1.y <= v2.y && v1.z <= v2.z)
		return true;
	else
		return false;
}
template<typename T>
bool operator>= (const Vector<T>& v1, const Vector<T>& v2)
{
	if(v1.x >= v2.x && v1.y >= v2.y && v1.z >= v2.z)
		return true;
	else
		return false;
}
template<typename T>
bool operator< (const Vector<T>& v1, const Vector<T>& v2)
{
	if(v1.x < v2.x)
		return true;
	if(v1.x == v2.x && v1.y < v2.y)
		return true;
	if(v1.x == v2.x && v1.y == v2.y && v1.z < v2.z)
		return true;
	return false;
}
template<typename T>
bool operator> (const Vector<T>& v1, const Vector<T>& v2)
{
	if(v1.x > v2.x && v1.y > v2.y && v1.z > v2.z)
		return true;
	else
		return false;
}

template<typename T>
T& Vector<T>::operator[](int index)
{
	switch(index)
	{
		case 0:
			return x;
			break;
		case 1:
			return y;
			break;
		case 2:
			return z;
			break;
		default:
			exit(0);
	}
}

#endif
