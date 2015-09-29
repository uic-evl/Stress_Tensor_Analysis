/* written by Abedul Haque, December 2009		*/
/* for CS3610 class project						*/

#pragma once

#include "emmintrin.h"

class TransferFunction{
private:
	float pointList[10][2] ;		// maximum of 10 points in each transfer function
	int index ; 
public:
	TransferFunction() ;
	void addPoint(float a, float b) ;
	float getValue(float a) ;
	__m128 getValue(__m128 m128_Intensity) ;
} ;