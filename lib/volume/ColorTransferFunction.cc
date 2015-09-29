/* written by Abedul Haque, December 2009		*/
/* for CS3610 class project						*/


#include"ColorTransferFunction.hh"

ColorTransferFunction::ColorTransferFunction(){
	index = -1 ;
} 

void ColorTransferFunction::addPoint(float a, float r, float g, float b){
	if( index >= 9 )
		return ;
	else{
		index++ ;
		this->pointList[index][0] = a ;
		this->pointList[index][1] = r ; 
		this->pointList[index][2] = g ; 
		this->pointList[index][3] = b ; 
	}
	return ;
} 


__m128 ColorTransferFunction::getValue(float a, float opacity){
	int i ;
	__m128 m128_pixelRGB ;
	for( i = 1 ; i < index ; i++){
		if ( a < pointList[i][0] )		// stop if you have hit a greater value
			break ;
	}

	int lowIndex = i-1 ;
	int highIndex = i ;

	// linear interpolate between these points	
	((float*)&m128_pixelRGB)[0]  = ( pointList[lowIndex][1] + (a-pointList[lowIndex][0])*
		 (pointList[highIndex][1] - pointList[lowIndex][1])/(pointList[highIndex][0] - pointList[lowIndex][0]) ) * opacity ;
	((float*)&m128_pixelRGB)[1] = ( pointList[lowIndex][2] + (a-pointList[lowIndex][0])*
		 (pointList[highIndex][2] - pointList[lowIndex][2])/(pointList[highIndex][0] - pointList[lowIndex][0]) ) * opacity ;
	((float*)&m128_pixelRGB)[2] = ( pointList[lowIndex][3] + (a-pointList[lowIndex][0])*
		 (pointList[highIndex][3] - pointList[lowIndex][3])/(pointList[highIndex][0] - pointList[lowIndex][0]) ) * opacity ;

	return m128_pixelRGB ;
} 

__m128 ColorTransferFunction::getColorValue(float a){
	int i ;
	__m128 m128_pixelRGB ;
	for( i = 1 ; i < index ; i++){
		if ( a < pointList[i][0] )		// stop if you have hit a greater value
			break ;
	}

	int lowIndex = i-1 ;
	int highIndex = i ;

	// linear interpolate between these points	
	((float*)&m128_pixelRGB)[0]  = ( pointList[lowIndex][1] + (a-pointList[lowIndex][0])*
		 (pointList[highIndex][1] - pointList[lowIndex][1])/(pointList[highIndex][0] - pointList[lowIndex][0]) ) ;
	((float*)&m128_pixelRGB)[1] = ( pointList[lowIndex][2] + (a-pointList[lowIndex][0])*
		 (pointList[highIndex][2] - pointList[lowIndex][2])/(pointList[highIndex][0] - pointList[lowIndex][0]) ) ;
	((float*)&m128_pixelRGB)[2] = ( pointList[lowIndex][3] + (a-pointList[lowIndex][0])*
		 (pointList[highIndex][3] - pointList[lowIndex][3])/(pointList[highIndex][0] - pointList[lowIndex][0]) ) ;

	return m128_pixelRGB ;
} 