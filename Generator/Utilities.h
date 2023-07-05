#ifndef UTIL_H
#define UTIL_H
#include "Body.h"

using namespace std;

class Util
{
	public:
	static Body make_mesh(Body, int);
	static Body merge_shapes(Body, Body);

	private:
	static Body make_mesh_handler(Body&, int, int);

};

#endif
