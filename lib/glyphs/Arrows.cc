#include "VectorFieldRenderer.hh"

/* vector field without eigen values */

void VectorFieldRenderer::Configure( Vector*** vf,iPoint volSize, Vector vSpacing ){
	
	this->VectorField = vf ;

	this->iSlices = volSize.z ;
	this->iRows = volSize.y ;
	this->iCols = volSize.x ;
	
	this->rendered = false;
}

void VectorFieldRenderer::DrawGlyph(){

	GLfloat angle = 0;
	GLfloat tip = 10.0 ;
	GLfloat rx, ry, prx, pry;
	GLfloat rz = 0.7*tip ;
	GLfloat radius = 1.0 ;

	GLfloat angleIncr = 45;

	glBegin(GL_LINES);
		glColor3f(1.0,1.0,1.0) ;
		glVertex3f(0,0,0) ;
		glVertex3f(0,0,tip) ;
	glEnd() ;
	// 
	// prx = radius*cos(MY_PI*angle/180) ;
	// pry = radius*sin(MY_PI*angle/180) ;
	// angle += angleIncr ;
	// glColor3f(0.9,0.6,0.6) ;
	// while( angle <= 360){
	// 	rx = radius*cos(MY_PI*angle/180) ;
	// 	ry = radius*sin(MY_PI*angle/180) ;		
	// 
	// 	glBegin(GL_LINE_LOOP) ;
	// 		glVertex3f(prx, pry, rz) ;
	// 		glVertex3f(rx, ry, rz) ;
	// 		glVertex3f(0,0,tip) ;
	// 	glEnd() ;
	// 
	// 	prx = rx ;
	// 	pry = ry ;
	// 	angle += angleIncr ;
	// }	
}
void VectorFieldRenderer::DrawAlignedGlyph(Vector v, Vector e){

	GLfloat angle1, angle2, temp ;
	glPushMatrix();
	
	//if( v.x != 0 )
	
	if( v.VectorZ == 0 && v.VectorX == 0 )
		angle1 = 0 ;
	else{
		
		angle1 = atan2(v.VectorX,v.VectorZ) ;
		angle1 = 180.0*angle1/MY_PI ;
		if( angle1 < 0 )
			angle1 += 360.0 ;
	}

	temp = sqrt(v.VectorZ*v.VectorZ + v.VectorX*v.VectorX) ;
	if( v.VectorY == 0 && temp == 0 )
		angle2 = 0 ;
	else{
		angle2 = atan2(v.VectorY,temp);
		angle2 = 180.0*angle2/MY_PI ;
		if (angle2 < 0.0 ){
			angle2 += 360.0 ;
		}
	}

	glRotatef(angle1,0,1,0) ;
	glRotatef(-1*angle2,1,0,0) ;

	DrawGlyph() ; 
	glPopMatrix() ;

}

void VectorFieldRenderer::setZoom(GLfloat zoom){
	
	GLfloat tempZoom = this->zoomScale; // preserve old value
	zoomScale = zoom; // set zoom
	RenderSlice(); // render grid
	zoomScale = tempZoom; // set scale back
}

void VectorFieldRenderer::RenderSlice(){
				
	int i , j , k , index;
	GLfloat posX, posY, posZ ;
	Vector dirV ;
	Vector eigV;

	cleanup();
	GLuint locInd;
	GLuint *locList;
	
	locInd = glGenLists(this->iRows*this->iCols*this->iSlices);
	locList = (GLuint *)malloc( sizeof(GLuint)*this->iRows*this->iCols*this->iSlices);
	
	int ind = 0;
	
	glColor3f(1.0, 1.0, 1.0); // white
	
	if( !this->bVisible )			// if not visible, not need to render
		return ;

	if( this->bSliceZVisible ) { 		

		k = this->iVectorFieldSliceNumberZ ;
		for(i = 0 ; i < this->iRows  ; i+=2){
			
			if( i%4 )
				continue ;
		
			glNewList(locInd+ind, GL_COMPILE);
			locList[ind] = ind;
			ind++;
			
			for ( j = 0 ; j < this->iCols ; j+=2){
				if (j%4)
					continue ;
				
				dirV = this->VectorField[k][i][j]  ;
				eigV = this->EigenValues[k][i][j] ;
				if( ( dirV.VectorX == dirV.VectorY) && ( dirV.VectorY == dirV.VectorZ ) && ( dirV.VectorZ == 0.0)  ){
						 posX = j*boxLenX + boxLenX/2 ;
						 posY = i*boxLenY + boxLenY/2 ;
						 posZ = k*boxLenZ + boxLenZ/2 ;
						 	
						 glPushMatrix() ;
						 glTranslatef(posX,posY,posZ) ;
						 glColor3f(0.4, 0.4, 0.4) ;
						 glutSolidSphere(2, 10, 10) ;
						 glPopMatrix() ; 
						continue ;				
					} 
				
				dirV.Normalize() ; 
				posX = j*boxLenX + boxLenX/2 ;
				posY = i*boxLenY + boxLenY/2 ;
				posZ = k*boxLenZ + boxLenZ/2 ;

				glPushMatrix() ;
				glTranslatef(posX,posY,posZ) ;
				glColor3f(1.0,1.0,1.0);
				glScalef(2.0,2.0,2.0);
				DrawSuperQuadric(eigV);
				//DrawAlignedGlyph(dirV, eigV) ;				
				glPopMatrix() ;
								
			} //for ( j = 0 ; j < this->iCols ; j++)
						
		} //for(i = 0 ; i < this->iRows  ; i++)
		
		glEndList();
		
	} //if( this->bSliceZVisible )	

	if( this->bSliceYVisible ){
		
		i = this->iVectorFieldSliceNumberY;
		for(k = 0 ; k < this->iSlices ; k+=2){
			if( k%4 )
				continue ;
				
			glNewList(locInd+ind, GL_COMPILE);
			locList[ind] = ind;
			ind++;
			
			
			for ( j = 0 ; j < this->iCols ; j+=2){
				if (j%4)
					continue ;
				index =  k*this->iRows*this->iCols + i * this->iCols + j ;

				dirV = this->VectorField[k][i][j] ; 
				eigV = this->EigenValues[k][i][j] ;
				
				 if( ( dirV.VectorX == dirV.VectorY) && ( dirV.VectorY == dirV.VectorZ ) && ( dirV.VectorZ == 0.0)  ){
					 posX = j*boxLenX + boxLenX/2 ;
								 posY = i*boxLenY + boxLenY/2 ;
								 posZ = k*boxLenZ + boxLenZ/2 ;
								 glPushMatrix() ;
								glTranslatef(posX,posY,posZ) ;
							 	glColor3f(0.4, 0.4, 0.4) ;
								glScalef(0.5,0.5,0.5);
								glutSolidSphere(2, 10, 10) ;
								glPopMatrix() ; 
								continue ;
							} 

				dirV.Normalize() ; 
				posX = j*boxLenX + boxLenX/2 ;
				posY = i*boxLenY + boxLenY/2 ;
				posZ = k*boxLenZ + boxLenZ/2 ;

				glPushMatrix() ;
				glTranslatef(posX,posY,posZ) ;
				glColor3f(1.0,1.0,1.0);
				glScalef(2.0,2.0,2.0);
				DrawSuperQuadric(eigV);				
				glPopMatrix() ;
			} //for ( j = 0 ; j < this->iCols ; j++)
						
		} //for(k = 0 ; k < this->iSlices ; k++)
	
		glEndList();
	
	} //else if( this->bSliceYVisible )
	if( this->bSliceXVisible ){
		//glPushMatrix();
		//glTranslatef(-0.5,0,0) ;
		j = this->iVectorFieldSliceNumberX  ;
		for(k = 0 ; k < this->iRows ; k+=2){
			
			if(k%4)
				continue ;
		
			glNewList(locInd+ind, GL_COMPILE);
			locList[ind] = ind;
			ind++;
			
			for ( i = 0 ; i < this->iRows ; i+=2){
				if (i%4)
					continue ;
				index =  k*this->iRows*this->iCols + i * this->iCols + j ;

				dirV = this->VectorField[k][i][j]  ; 
				eigV = this->EigenValues[k][i][j] ;
				
				if( ( dirV.VectorX == dirV.VectorY) && ( dirV.VectorY == dirV.VectorZ ) && ( dirV.VectorZ == 0.0)  ){
										
					posX = j*boxLenX + boxLenX/2 ;
					posY = i*boxLenY + boxLenY/2 ;
					posZ = k*boxLenZ + boxLenZ/2 ;
					glPushMatrix() ;
					glTranslatef(posX,posY,posZ) ;
					glColor3f(0.2, 0.2, 0.2) ;
					glScalef(0.5,0.5,0.5);
					glutSolidSphere(2, 10, 10) ;
					glPopMatrix() ; 
					continue ;
				} 

				dirV.Normalize() ; 
				posX = j*boxLenX + boxLenX/2 ;
				posY = i*boxLenY + boxLenY/2 ;
				posZ = k*boxLenZ + boxLenZ/2 ;

				glPushMatrix() ;
				glTranslatef(posX,posY,posZ) ;
				glColor3f(1.0,1.0,1.0);
				glScalef(2.0,2.0,2.0);
				//DrawSuperQuadric(eigV);
				DrawAlignedGlyph(dirV, eigV) ;				
				glPopMatrix() ;
			} //for ( i = 0 ; i < this->iRows ; i++)
			
			glEndList();
			
		} //for(k = 0 ; k < this->iSlices ; k++)		
	
	} //else if( this->bSliceXVisible )
		
	this->lists = locList;
	this->index = locInd;	
	
	this->rendered = true;
}

void VectorFieldRenderer::IncreaseSliceNumber(){
		
	if( this->bSliceXVisible  ){
		this->iVectorFieldSliceNumberX++ ;
		this->iVectorFieldSliceNumberX = min(this->iVectorFieldSliceNumberX, this->iCols-1 ) ;
	}		  
	if( this->bSliceYVisible ){
		this->iVectorFieldSliceNumberY++ ;
		this->iVectorFieldSliceNumberY = min(this->iVectorFieldSliceNumberY, this->iRows-1 ) ;
	}
	if( this->bSliceZVisible ){
		this->iVectorFieldSliceNumberZ++ ;
		this->iVectorFieldSliceNumberZ = min(this->iVectorFieldSliceNumberZ, this->iSlices-1 ) ;
	}
} 

void VectorFieldRenderer::DecreaseSliceNumber(){
		
	if( this->bSliceXVisible  ){
		this->iVectorFieldSliceNumberX-- ;
		this->iVectorFieldSliceNumberX = max(this->iVectorFieldSliceNumberX, 0) ;
	}		  
	if( this->bSliceYVisible ){
		this->iVectorFieldSliceNumberY-- ;
		this->iVectorFieldSliceNumberY = max(this->iVectorFieldSliceNumberY, 0) ;
	}
	if( this->bSliceZVisible ){
		this->iVectorFieldSliceNumberZ-- ;
		this->iVectorFieldSliceNumberZ = max(this->iVectorFieldSliceNumberZ, 0) ;
	}
} 

void VectorFieldRenderer::ToggleSliceVisibilityX(){
	if( this->bSliceXVisible  )
		this->bSliceXVisible   = false ;		  
	else{
		this->bSliceXVisible  = true ;		  
	}
} 

void VectorFieldRenderer::ToggleSliceVisibilityY(){
	if( this->bSliceYVisible  )
		this->bSliceYVisible   = false ;		  
	else{
		this->bSliceYVisible  = true ;		  
	}
} 

void VectorFieldRenderer::ToggleSliceVisibilityZ(){
	if( this->bSliceZVisible  )
		this->bSliceZVisible   = false ;		  
	else{
		this->bSliceZVisible  = true ;		  
	}
}