#ifndef VECTOR3_INCLUDED
#define VECTOR3_INCLUDED

#include "Point3.hh"

class Vector3 {

public:	
	float x, y, z;

	Vector3(){ x = y = z = 0 ; }
	Vector3(float xx, float yy, float zz){ x = xx ; y = yy ; z = zz ;}
	Vector3( Vector3 &v){ x = v.x ; y = v.y ; z = v.z ; }

	void set(float dx, float dy, float dz){ x = dx; y = dy; z = dz ;}	
	void set(Vector3 &v){ x = v.x ; y = v.y ; z = v.z ; }

	void flip(){ x = -x ; y = -y ; z = -z ; }
	void setDiff( Point3 &a, Point3 &b){ x = a.x  - b.x ; y = a.y - b.y ; z = a.z - b.z ; }

	void normalize() ;
	Vector3 cross(Vector3 b) ;	
	float dot(Vector3 b );

};

#endif