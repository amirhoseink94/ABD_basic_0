#ifndef UTIL_H
#define UTIL_H
#include "Shape.h"

using namespace std;

class Util
{
	public:
	static Shape make_mesh(Shape, int);
	static Shape merge_shapes(Shape, Shape);

	private:
	static Shape make_mesh_handler(Shape, int, int);
	
};

#endif
