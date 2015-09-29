/* written by Abedul Haque, December 2009		*/
/* for CS3610 class project						*/

#include "TransferFunction.hh"

TransferFunction::TransferFunction(){
	index = -1 ;
} 

void TransferFunction::addPoint(float a, float b){
	if( index >= 9 )
		return ;
	else{
		index++ ;
		this->pointList[index][0] = a ;
		this->pointList[index][1] = b ; 
	}
	return ;
} 

float TransferFunction::getValue(float a){
	int i ;

	for( i = 1 ; i < index ; i++){
		if ( a < pointList[i][0] )		// stop if you have hit a greater value
			break ;
	}

	int lowIndex = i-1 ;
	int highIndex = i ;

	// linear interpolate between these points	
	return		pointList[lowIndex][1] +
				(a-pointList[lowIndex][0])*
				(pointList[highIndex][1] - pointList[lowIndex][1])/(pointList[highIndex][0] - pointList[lowIndex][0]) ;
} 

__m128 TransferFunction::getValue(__m128 m128_Intensity){

	int point ;
	__m128 m128_Opacity, m128_LowValX, m128_HighValX, m128_LowValY, m128_HighValY ;
	
	for(int i = 0 ; i < 4 ; i++){			// process for each component	
		float a = (float)( ((float*)&m128_Intensity)[i] ) ; 
		for( point = 1 ; point < index ; point++){		// find low and high anchor point for interpolation
			if ( a <= pointList[point][0] )		// stop if you have hit a greater value
				break ;
		}

		((float*)&m128_LowValX)[i] =  pointList[point-1][0] ;
		((float*)&m128_HighValX)[i] =  pointList[point][0] ;

		((float*)&m128_LowValY)[i] =  pointList[point-1][1] ;
		((float*)&m128_HighValY)[i] =  pointList[point][1] ;
	}

	m128_Opacity = _mm_add_ps( m128_LowValY,
						_mm_mul_ps(
						_mm_sub_ps(m128_Intensity, m128_LowValX),
						_mm_div_ps( _mm_sub_ps(m128_HighValY, m128_LowValY),
									_mm_sub_ps(m128_HighValX, m128_LowValX) )
						)
					);
	return m128_Opacity ;
} 
