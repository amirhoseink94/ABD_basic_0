#include "../include/Vector.h"
#include <cmath>

template<typename T>
Vector<T>::Vector()
{
	this->x = 0;
	this->y = 0;
	this->z = 0;
}
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
    temp.x = this->x + obj.x;
    temp.y = this->y + obj.y;
    temp.z = this->z + obj.z;
    return temp;
}
template<typename T>
Vector<T> Vector<T>::operator-(const Vector<T>& obj)
{
    Vector<T> temp;
    temp.x = this->x - obj.x;
    temp.y = this->y - obj.y;
    temp.z = this->z - obj.z;
    return temp;
}
template<typename T>
Vector<T>& Vector<T>::operator= (const Vector<T>& obj)
{
	this->x = obj.x;
	this->y = obj.y;
	this->z = obj.z;
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
	if(v1.x < v2.x && v1.y < v2.y && v1.z < v2.z)
		return true;
	else
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
