#ifndef POINT3_INCLUDED
#define POINT3_INCLUDED

class Point3{
public:
	float x, y, z;
	Point3(){ x = y = z = 0.0f ;}
	Point3(float xx, float yy, float zz ){ x = xx ; y = yy ; z = zz ;}

	void set( float xx, float yy, float zz ){ x = xx; y = yy ; z = zz ;}
	void set( Point3 &p ){ x = p.x ; y = p.y ; z = p.z ;}

	//void build4tuple()
} ;

/*class PointStack{
public:
	PointStack(){ stack = NULL ; } 
	void push( Point3 p) ;
	Point3 pop() ;
	void release() ;

private :
	typedef struct pStack{
		Point3 *current ;
		struct pStack *next ;
	};
	pStack *stack ;
}; */

#endif
