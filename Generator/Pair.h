#ifndef PAIR_H
#define PAIR_H

#include <iostream>
#include <cmath>
using namespace std;

template<typename T>
class Pair
{
	public:
    T x, y;
    Pair();
    Pair(T x, T y);
    // Overload the * operator
    Pair<T> operator+ (const Pair<T>& obj);
    Pair<T> operator- (const Pair<T>& obj);
    Pair<T>& operator= (const Pair<T>& obj);
    T& operator[](int);
		//const T operator[](int);
    float length();
    Pair<T> num_multi(float);

    template <typename U>
    friend ostream& operator<<(ostream& os, const Pair<U>& obj);
    template <typename U>
    friend bool operator== (const Pair<U>& v1, const Pair<U>& v2);
    template <typename U>
    friend bool operator!= (const Pair<U>& v1, const Pair<U>& v2);
    template <typename U>
    friend bool operator<= (const Pair<U>& v1, const Pair<U>& v2);
    template <typename U>
    friend bool operator>= (const Pair<U>& v1, const Pair<U>& v2);
    template <typename U>
    friend bool operator< (const Pair<U>& v1, const Pair<U>& v2);
    template <typename U>
    friend bool operator> (const Pair<U>& v1, const Pair<U>& v2);

};



template<typename T>
Pair<T>::Pair()
{

}
template<typename T>
Pair<T>::Pair(T x, T y)
{
	this->x = x;
	this->y = y;
}
template<typename T>
Pair<T> Pair<T>::operator+(const Pair<T>& obj)
{
    Pair<T> temp;
    temp.x = x + obj.x;
    temp.y = y + obj.y;
    return temp;
}
template<typename T>
Pair<T> Pair<T>::operator-(const Pair<T>& obj)
{
    Pair<T> temp;
    temp.x = x - obj.x;
    temp.y = y - obj.y;
    return temp;
}
template<typename T>
Pair<T>& Pair<T>::operator= (const Pair<T>& obj)
{
	x = obj.x;
	y = obj.y;
    return *this;
}
template<typename T>
ostream& operator<<(ostream& os, const Pair<T>& obj)
{
    os << "(" << obj.x << ", " << obj.y  << ")";;
    return os;
}
template<typename T>
Pair<T> Pair<T>::num_multi(float a)
{
	Pair<T> temp;
    temp.x = a * x;
    temp.y = a * y;

    return temp;
}
template<typename T>
float Pair<T>::length()
{
	float l = pow(x, 2) + pow(y, 2);
	l = sqrt(l);
	return l;
}
template<typename T>
bool operator== (const Pair<T>& v1, const Pair<T>& v2)
{
	if(v1.x == v2.x && v1.y == v2.y)
		return true;
	else
		return false;
}
template<typename T>
bool operator!= (const Pair<T>& v1, const Pair<T>& v2)
{
	if(v1.x != v2.x || v1.y != v2.y)
		return true;
	else
		return false;
}
template<typename T>
bool operator<= (const Pair<T>& v1, const Pair<T>& v2)
{
	if(v1.x <= v2.x && v1.y <= v2.y)
		return true;
	else
		return false;
}
template<typename T>
bool operator>= (const Pair<T>& v1, const Pair<T>& v2)
{
	if(v1.x >= v2.x && v1.y >= v2.y)
		return true;
	else
		return false;
}
template<typename T>
bool operator< (const Pair<T>& v1, const Pair<T>& v2)
{
	if(v1.x < v2.x)
		return true;
	if(v1.x == v2.x && v1.y < v2.y)
		return true;
	return false;
}
template<typename T>
bool operator> (const Pair<T>& v1, const Pair<T>& v2)
{
	if(v1.x > v2.x && v1.y > v2.y)
		return true;
	else
		return false;
}

template<typename T>
T& Pair<T>::operator[](int index)
{
	switch(index)
	{
		case 0:
			return x;
			break;
		case 1:
			return y;
			break;
		default:
			exit(0);
	}
}
/*template<typename T>
const T Pair<T>::operator[](int index)
{
	switch(index)
	{
		case 0:
			return x;
			break;
		case 1:
			return y;
			break;
		default:
			exit(0);
	}
}*/

#endif
