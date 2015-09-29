/* written by Abedul Haque, December 2009		*/
/* for CS3610 class project						*/


#pragma once
#include "emmintrin.h"


class ColorTransferFunction{
private:
	float pointList[10][4] ;		// maximum of 10 points in each transfer function
	int index ; 
public:
	ColorTransferFunction() ;
	void addPoint(float a, float r, float g, float b) ;
	__m128 getValue(float a, float opacity) ;
	__m128 getColorValue(float a) ;
} ;