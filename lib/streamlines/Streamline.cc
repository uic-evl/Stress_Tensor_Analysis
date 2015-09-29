#include"Streamline.hh"

extern GLuint createVBO(GLfloat *points, int size);
extern GLuint compileShader(const char *vsource, const char *fsource, const char *gsource);
extern void deleteVBO(GLuint vbo);
extern void draw_streamlines(GLuint vbo, GLuint shader, int size, GLfloat color[3]);
extern char *textFileRead(char *fn) ;

Streamline::Streamline(){
	this->bReadyToProcess = false ;		// We are not ready to generate streamlines
	this->bVisible = false	;			// streamlines are not visible			
	TempSelectedStreamline = -1;		// no streamline temporarily selected

	//three selected streamline color
	streamlineColor[0][0] = 0.2 ; streamlineColor[0][1] = 0.3 ;	streamlineColor[0][2] = 0.8 ;	
	streamlineColor[1][0] = 0.3 ; streamlineColor[1][1] = 0.8 ;	streamlineColor[1][2] = 0.4 ;	
	streamlineColor[2][0] = 0.8 ; streamlineColor[2][1] = 0.3 ;	streamlineColor[2][2] = 0.2 ;	
	//unselected streamline
	streamlineColor[3][0] = 0.4 ; streamlineColor[3][1] = 0.4 ;	streamlineColor[3][2] = 0.4 ;	
	//temporary selected streamline
	streamlineColor[4][0] = 0.5 ; streamlineColor[4][1] = 0.5 ;	streamlineColor[4][2] = 0.5 ;		
} 

void Streamline::SetSpacing(float dx, float dy, float dz){
	this->boxLenX = dx ;
	this->boxLenY = dy ;
	this->boxLenZ = dz ;
}

void Streamline::biLinearInterpolate(Vector &v1, Vector &v2, Vector &v3, Vector &v4, float fDelta1, float fDelta2){
	v1.Interpolate(v2,fDelta1) ;
	v3.Interpolate(v4,fDelta1) ;
	v1.Interpolate(v3,fDelta2) ; 
}

void Streamline::Configure(fPoint SeedPoint, int n){
	this->fSeedPosX = SeedPoint.x ;								// set the seed point
	this->fSeedPosY = SeedPoint.y ;
	this->fSeedPosZ = SeedPoint.z ;
	this->iNumberOfStreamlines = min(MAX_NUMBER_STREAMLINE, n) ;//# of streamlines needed

	StreamlineDisplayListIndex = glGenLists(this->iNumberOfStreamlines) ;		//Creating display list indices
	if( StreamlineDisplayListIndex == 0 ){
		printf("Can't allocate display list indices.\n") ;
		this->bReadyToProcess = false ;
	}
	else
		this->bReadyToProcess = true ;							//All data has been set. We can generate streamline now

	for( int i = 0 ; i < this->iNumberOfStreamlines ; i++){
		SelectedStreamlines[i] = UNSELECTED ;					// not selected
		streamlineColorIndex[i] = 3 ;							// not selected color
	}
	for(int i = 0 ; i < 3 ; i++ )								
		this->freeColorForSelected[i] = 1 ;						// this color has not been used		

	this->bSeedpointSelected = false ;
	this->bSeedpointMoved = false ;
}

void Streamline::DrawSeedPoint(){
	glPushMatrix() ;
	glTranslatef(this->fSeedPosX, this->fSeedPosY, this->fSeedPosZ) ;   
	if( this->bSeedpointSelected )
		glColor3f(0.5, 0.5, 0.0) ;
	else
		glColor3f(0.0, 0.5, 0.5) ;
	glutSolidSphere(3, 10, 10) ;
	glPopMatrix() ;
} 

// once you delete strealines you have to configure again with all the data
void Streamline::DeleteStreamlines(){
	glDeleteLists(this->StreamlineDisplayListIndex, this->iNumberOfStreamlines) ;	// delete display list
	this->bReadyToProcess = false ;													// must reconfigure to make it true
} 

// returns false if the position is out of bound
bool Streamline::DerivativeAtPosition(Vector CurrPos, Vector& vDerivative){
	int iXMin = floor(CurrPos.X()) ;
	int iXMax = ceil (CurrPos.X()) ;		
	float fDelta1 = CurrPos.X() - (float)iXMin ;

	int iYMin = floor(CurrPos.Y()) ;
	int iYMax = ceil(CurrPos.Y()) ;
	float fDelta2 = CurrPos.Y() - (float)iYMin;

	int iZMin = floor(CurrPos.Z()) ;
	int iZMax = ceil(CurrPos.Z()) ;
	float fDelta3 = CurrPos.Z() - (float)iZMin ;

	if( iXMin < 0 || iXMin >= this->cols  )			// condition check. if your current position is outside
		return false;									// the volume, then don't continue with the line
	if( iXMax < 0 || iXMax >= this->cols )				
		return false;		
	if( iYMin < 0 || iYMax >= this->rows  )			// I am not sure we need all the conditions
		return false;									// may be we can do Min with 0 and Max with ROWS/COLS/SLICES
	if( iYMax < 0 || iYMax >= this->rows )				// but right now I am not taking any chance
		return false;	
	if( iZMin < 0 || iZMin >= this->slices )
		return false;	
	if( iZMax < 0 || iZMax >= this->slices )
		return false;	

	// get four corners of the first XY plane
	//int index = iZMin*ROWS*COLS + iYMin*rows +
	Vector v1 =  vectorField[iZMin][iYMin][iXMin] ;
	Vector v2 =  vectorField[iZMin][iYMin][iXMax] ;
	Vector v3 =  vectorField[iZMin][iYMax][iXMin] ;
	Vector v4 =  vectorField[iZMin][iYMax][iXMax] ;
	biLinearInterpolate(v1, v2, v3, v4, fDelta1, fDelta2) ;		// doing bi linear interpolation in this plane	
																				// result is in v1
	// get four corners of the first XY plane
	//int index = iZMin*ROWS*COLS + iYMin*rows +
	Vector v5 =  vectorField[iZMax][iYMin][iXMin] ;
	Vector v6 =  vectorField[iZMax][iYMin][iXMax] ;
	Vector v7 =  vectorField[iZMax][iYMax][iXMin] ;
	Vector v8 =  vectorField[iZMax][iYMax][iXMax] ;
	biLinearInterpolate(v5, v6, v7, v8, fDelta1, fDelta2) ;		// doing bi linear interpolation in this plane	
																				// result is in v5	
					
	// now interpolate v1 and v5
	v1.Interpolate(v5, fDelta3) ;
	v1.Normalize() ;
	vDerivative.VectorX = v1.X() ;
	vDerivative.VectorY = v1.Y() ;
	vDerivative.VectorZ = v1.Z() ;

	return true ;			// return success
} 
// returns false if the position is out of bound
// true if the position is good. Then new postion is set in NextPos vector
bool Streamline::EulerStep(Vector CurrPos, Vector& NextPos, float fStepSize){

	Vector vDerivative(0.0f, 0.0f, 0.0f) ;

	if( !this->DerivativeAtPosition(CurrPos, vDerivative) )			// if the position is out of bound
		return false ;

	vDerivative.ScalarMult(fStepSize) ; 
	NextPos.VectorX = CurrPos.X() + vDerivative.X() ;
	NextPos.VectorY = CurrPos.Y() + vDerivative.Y() ;
	NextPos.VectorZ = CurrPos.Z() + vDerivative.Z() ;

	return true;
}

bool Streamline::MidPointMethod(Vector CurrPos, Vector& NextPos, float fStepSize){

	Vector vDerivative(0.0f, 0.0f, 0.0f) ;

	if( !this->DerivativeAtPosition(CurrPos, vDerivative) )			// if the position is out of bound
		return false ;

	// now get the midpoint
	// k1 = h*Derivative
	// we are doing fStepSize/2 to get k1/2
	// Midpoint will have X0 + k1/2 which is mid point between the current position 
	// and next position according to euler step
	vDerivative.ScalarMult(fStepSize/2 ); 
	Vector MidPoint(0.0, 0.0, 0.0) ;
	MidPoint.VectorX = CurrPos.X() + vDerivative.X() ;
	MidPoint.VectorY = CurrPos.Y() + vDerivative.Y() ;
	MidPoint.VectorZ = CurrPos.Z() + vDerivative.Z() ;

	// Now get the derivative at midPoint
	if( !this->DerivativeAtPosition(MidPoint, vDerivative) )			// if the position is out of bound
		return false ;
	// Now find the next position according to this derivative
	vDerivative.ScalarMult(fStepSize);									// now take the full step
	NextPos.VectorX = CurrPos.X() + vDerivative.X() ;
	NextPos.VectorY = CurrPos.Y() + vDerivative.Y() ;
	NextPos.VectorZ = CurrPos.Z() + vDerivative.Z() ;

	return true;
}
bool Streamline::RK4Method(Vector CurrPos, Vector& NextPos, float fStepSize){

	Vector vK1, vK2, vK3, vK4;
	
	// finding K1
	if( !this->DerivativeAtPosition(CurrPos, vK1) )					// if the position is out of bound
		return false ;	
	vK1.ScalarMult(fStepSize) ;										// k1 = h*Derivative
	// now prepare for K2
	Vector vTemp1 = vK1 ;											// (k1)/2
	vTemp1.ScalarMult(0.5) ; 	
	Vector HalfK1Point;

	HalfK1Point.VectorX = CurrPos.X() + vTemp1.X() ;		// HalfK1Point = X0 + K1/2
	HalfK1Point.VectorY = CurrPos.Y() + vTemp1.Y() ;
	HalfK1Point.VectorZ = CurrPos.Z() + vTemp1.Z() ;
	// Now get the derivative at HalfK1Point which will be K2
	if( !this->DerivativeAtPosition(HalfK1Point, vK2) )				// calculating derivative at HalfK1 i.e f(X0+K1/2)
		return false ;
	
	vK2.ScalarMult(fStepSize) ;										// k2 = h*Derivative

	// Now prepare for K3
	Vector vTemp2 = vK2 ;											// vTemp = (k2)/2
	vTemp2.ScalarMult(0.5) ; 	
	Vector HalfK2Point;											
	HalfK2Point.VectorX = CurrPos.X() + vTemp2.X() ;		// HalfK2Point = X0 + K2/2
	HalfK2Point.VectorY = CurrPos.Y() + vTemp2.Y() ;
	HalfK2Point.VectorZ = CurrPos.Z() + vTemp2.Z() ;
	// Now get the derivative at HalfK2Point
	if( !this->DerivativeAtPosition(HalfK2Point, vK3) )				// calculating derivative i.e. f(X0+K2/2)
		return false ;

	vK3.ScalarMult(fStepSize) ;										// k3 = h*Derivative

	// Now prepare for K4
	Vector vTemp3 = vK3 ;											// vTemp = (K3), Note, there is no dividing by 2 here
	Vector K3Point;											
	K3Point.VectorX = CurrPos.X() + vTemp3.X() ;			// K3Point = X0 + K3
	K3Point.VectorY = CurrPos.Y() + vTemp3.Y() ;
	K3Point.VectorZ = CurrPos.Z() + vTemp3.Z() ;
	// Now get the derivative at K3Point
	if( !this->DerivativeAtPosition(K3Point, vK4) )					// calculating derivative i.e. f(X0+K3) 
		return false ;
	vK4.ScalarMult(fStepSize) ;										// k4 = h*Derivative

	// Now find the next position according to all these derivates
	vK1.ScalarMult((1.0/6.0)) ;
	vK2.ScalarMult((1.0/3.0));
	vK3.ScalarMult((1.0/3.0)) ;
	vK4.ScalarMult((1.0/6.0)) ;

	NextPos.VectorX = CurrPos.X() ;
	NextPos.VectorY = CurrPos.Y() ;
	NextPos.VectorZ = CurrPos.Z() ;

	NextPos.Add(vK1) ;						// NextPos = CurrPos + (1/6)K1
	NextPos.Add(vK2);						// NextPos = CurrPos + (1/6)K1 + (1/3)K2
	NextPos.Add(vK3);						// NextPos = CurrPos + (1/6)K1 + (1/3)K2 + (1/3)K3 
	NextPos.Add(vK4);						// NextPos = CurrPos + (1/6)K1 + (1/3)K2 + (1/3)K3 + (1/6)K4

	return true;
}
void Streamline::CreateStreamlines(){		
	char *g = textFileRead("/Users/onorinbejasus/Dropbox/Vis/Projects/Tensor/lib/draw/geometry.geom");
	shader = compileShader(vertexShader, linePixelShader, g);

	Vector CurrPos, NextPos ;

	int iStepCount = 0 ;					// integration time step
	float fStepSize = 0.25 ;				// stepsize of each step
	int iMaxStepCount = 0 ;					// just to record the maximum number of steps required

	int ii  = this->StreamlineDisplayListIndex ; // display list index for the first streamline
	
	int index = 0;

	for( int k = -2 ; k <= 2 ; k+=2){
		for(int i = -2 ; i <= 2 ; i+=2){
			for(int j = -2 ; j <= 2 ; j+=2){
				glNewList(ii, GL_COMPILE);
				CurrPos.VectorX = fSeedPosX + j ;
				CurrPos.VectorY = fSeedPosY + i ;
				CurrPos.VectorZ = fSeedPosZ + k;
				iStepCount = 0 ;

				while( iStepCount < 40000 )			// max of 10000 steps
				{
					
					// do one RK4 step
					if( !this->RK4Method (CurrPos, NextPos, fStepSize) )		
						break ;
					
					stream_pts.push_back(new Vector(CurrPos.X()*boxLenX, 
								CurrPos.Y()*boxLenY,CurrPos.Z()*boxLenZ));
					
					points.push_back( CurrPos.X());
					points.push_back( CurrPos.Y());
					points.push_back( CurrPos.Z());
					points.push_back(1.0);
					
					glBegin(GL_LINES) ;											// Now, draw the line
						glVertex3f(CurrPos.VectorX*this->boxLenX, CurrPos.VectorY*this->boxLenY, CurrPos.VectorZ*this->boxLenZ) ;
						glVertex3f(NextPos.VectorX*this->boxLenX, NextPos.VectorY*this->boxLenY, NextPos.VectorZ*this->boxLenZ) ;
					glEnd() ;
					
					CurrPos = NextPos ;		// update current					
					iStepCount++ ;				//update stepCount
				} //while( iStepCount < 10000 )
				
				glEndList() ;
				
				// create the vbo				
				vbo[index] = createVBO(&points[0], points.size());
				sizes[index++] = points.size();

				points.clear();

				ii++ ;
				if( iMaxStepCount < iStepCount )
					iMaxStepCount = iStepCount ;
			} //for(int j = -2 ; j <= 2 ; j+=2)
		} //for(int i = -2 ; i <= 2 ; i+=2)
	} //for( int k = -2 ; k <= 2 ; k+=2)

	//printf("Max steps required %d\n",iMaxStepCount) ;
} 

void Streamline::RenderStreamlines(GLfloat *eigenfloats){
	int iDisplayListIndex ;
	
	if( this->bVisible ){								// if they are visible then draw. Otherwise no need.
		this->DrawSeedPoint() ; 
		int temp = 4;
		for(int i = 0 ; i < this->iNumberOfStreamlines ; i++){
			if(i % 3 == 0)
				temp = 4;//i % 4;
			else
				temp = 4;
			draw_streamlines(vbo[i], shader, sizes[i] * 4.0 * sizeof(GLfloat), streamlineColor[temp]);
		}
	}
} 

void Streamline::SelectModeSeedPointDraw(){
	if( !this->bVisible  )
		return ;
	glPushMatrix() ;
	glPushName(1) ;				// object name for seed point
	glTranslatef(this->fSeedPosX, this->fSeedPosY, this->fSeedPosZ) ;   
	glutSolidSphere(3, 10, 10) ;
	glPopName() ;
	glPopMatrix() ;
} 

void Streamline::SelectModeStreamlinesDraw(){

	int iDisplayListIndex ;
	if( !this->bVisible  )
		return ;
	glPushName(-1) ;					// placeholder for the next drawings
	for(int i = 0 ; i < this->iNumberOfStreamlines ; i++){
		iDisplayListIndex = this->StreamlineDisplayListIndex + i ; 
		glLoadName(iDisplayListIndex) ;					// display list index will work as a name for the object
		// printf("select mode display list %d\n",iDisplayListIndex) ;
		glCallList(iDisplayListIndex) ;
	}
	glPopName() ;
} 

void Streamline::UpdateStreamlineStatus(int iName){
	int iArrayIndex = iName - this->StreamlineDisplayListIndex ;
	if( this->SelectedStreamlines[iArrayIndex]  ){					// toggle status
		this->SelectedStreamlines[iArrayIndex] = UNSELECTED ;		// make it unselected
		int colorIndex = this->streamlineColorIndex[iArrayIndex] ;	// change the color	
		this->freeColorForSelected[colorIndex] = 1 ;				// make the color free 
		this->streamlineColorIndex[iArrayIndex] =  3 ;		// change the color	of the streamline
	}
	else{
		this->SelectedStreamlines[iArrayIndex] = SELECTED ;		
		//find a free color
		int freeColor = 0 ;
		for( int i = 0 ; i < 3 ; i++){
			if( this->freeColorForSelected[i] ){
				freeColor = i ;
				break ;
			}
		}
		this->streamlineColorIndex[iArrayIndex] =  freeColor ;
		this->freeColorForSelected[freeColor] = 0 ;						// this color is no longer free 
	}
} 

void Streamline::TemporarySelectedStreamline( int iName){
	
	printf("streamline selected!\n");
	
	int iArrayIndex = iName - this->StreamlineDisplayListIndex ;
	this->TempSelectedStreamline = iArrayIndex ;
} 

void Streamline::UpdateSeedPoint(){
	this->DeleteStreamlines() ; 
	this->Configure( fPoint(this->fSeedPosX, this->fSeedPosY, this->fSeedPosZ) ,this->iNumberOfStreamlines) ;
	this->CreateStreamlines() ; 
} 

void Streamline::GetSeedPoint(fPoint& fp){
	//fPoint seed(this->fSeedPosX, this->fSeedPosY, this->fSeedPosZ ) ;
	fp.x = this->fSeedPosX ;
	fp.y = this->fSeedPosY ;
	fp.z = this->fSeedPosZ ;	
}
 
// change the seed point location but don't update streamlines 
void Streamline::MoveSeedPoint(fPoint& fp){
	this->fSeedPosX = fp.x ;
	this->fSeedPosY = fp.y ;
	this->fSeedPosZ = fp.z ;
	this->bSeedpointMoved = true ;
	printf("%f %f %f\n",this->fSeedPosX, this->fSeedPosY, this->fSeedPosZ ) ;
} 