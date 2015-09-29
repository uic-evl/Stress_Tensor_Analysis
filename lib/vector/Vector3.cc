#include <math.h>
#include "Vector3.hh"

void Vector3::normalize(){
	float length ;
	length = sqrt(x*x + y*y + z*z ) ;

	x = x/length ;
	y = y/length ;
	z = z/length ;
}

Vector3 Vector3::cross( Vector3 b ) {
	Vector3 normVector ;

	normVector.set( (y*b.z-z*b.y), (b.x*z-x*b.z), ( x*b.y-y*b.x) ) ;
	return normVector ;
}

float Vector3::dot( Vector3 b ){
	float temp ;
	temp = ( x*b.x + y*b.y + z*b.z );
	return temp ;
}




