/* written by Abedul Haque, December 2009		*/
/* for CS3610 class project						*/

#include "VolumeRenderer.hh"

void VolumeRenderer::Initialize(std::string szp, int width, int type){ 

	this->params.LoadParameters(szp);							// load the parameter object   
	this->volumeImg.LoadCtSansInterpolation(params.CtFileName,				// volume filename   
                                   params.CtFileSize,						// volume filesize
                                   params.CtFileVoxelSpacing,				// voxel size
                                   	params.CtSamplingInterval);				// sampling interval

   	volumeImg.vPosition    = new Vector(0, 0, 0);				// translation and
	volumeImg.vOrientation = new Vector(0, 0, 0);				// rotation of the volume
	volumeImg.ComputeTransforms();								// calculate the Vol-to-world and world-to-Vol transform matrices
	this->projectionWidth = width;
	this->projectionHeight = width ;
	this->InitializeProjector(type) ;											// set projector specs
	this->BuildColorBar() ; 

}

void VolumeRenderer::CreateVolImage()
{		
	this->volumeImg.ComputeTransforms();						// calculate the Vol-to-world and world-to-Vol transform matrices	
	this->projector.CalculateTransformations() ;				// update projector to world transformation matrices
	this->projector.Clear() ;									// clear the projector image
	if( this->bIntensity ){
		this->projector.RayCast(volumeImg);
		this->projector.Normalize() ; 
	}
	else{
		this->projector.ColorRayCast(volumeImg);
		this->projector.NormalizeRGB() ;
	}	
}

void VolumeRenderer::InitializeProjector(int type){

	Vector projectorLocation(0,0,0);

		projectorLocation = new Vector(0,0,2500) ;

	Vector projectorDirection(0.0, 0.0, -1.0) ;
	Vector centralSpot((double)projectionWidth/2,(double)projectionHeight/2, 0.0) ;
	iPoint viewportSize(projectionWidth,projectionHeight,0) ;
	float fSourceToDetector = this->params.RedSourceToDetector ;
	float fDetectorPixelSize = 0.125;
    this->projector.Configure(projectorLocation,                  // set up the projector system
                    projectorDirection,
                    fSourceToDetector,
                    viewportSize,
                    fDetectorPixelSize,
                    centralSpot);
    this->projector.iFastPassSamplingInterval = 8;				// configure the fast-pass projection method		
																// every 8th pixel is ray-casted. Rest are interpolated
	//updating the transfer functions
	this->projector.opacityFunction = this->params.opTranFunc ;
	this->projector.colorTranFun = this->params.colTranFunc ; 

	// the scaling factors
	this->projector.vProjectorScale.fVector[0] = 1;  
	this->projector.vProjectorScale.fVector[1] = 1;  
	this->projector.vProjectorScale.fVector[2] = 1;  
	//calculate transformation matrices
	this->projector.CalculateTransformations() ; 
	//projector.printProjector() ; 
}

void VolumeRenderer::RotateObject(float fAngleX, float fAngleY, float fAngleZ){
	this->bRayCastProjChanged = true ;			// we need to do new projection to get updated picture
	volumeImg.vOrientation = new Vector(fAngleX, fAngleY, fAngleZ) ; 
}  

void VolumeRenderer::ChangeProjectorLocation(float fDelta){
	this->bRayCastProjChanged = true ;			// we need to do new projection to get updated picture
	this->projector.vSourceLocation.VectorZ += fDelta  ; 
} 

void VolumeRenderer::RenderVolImage(){ 	
	if( this->bVisible){													// if volume image was enabled
		if( this->bRayCastProjChanged ){
			this->CreateVolImage() ;										// update the image 
			this->bRayCastProjChanged = false ;
		}

	    glRasterPos3f(	-1.0 * (float)projectionWidth/2   * this->projector.fDetectorPixelSize ,			// this is basically the position of the near clipping
						-1.0 * (float)projectionHeight/2  * this->projector.fDetectorPixelSize ,				// plane in object coordinate system	
						this->projector.vSourceLocation.VectorZ - this->projector.fSourceToDetector);
		if( this->bIntensity )
			glDrawPixels(this->projector.ipImageSize.y, 
						this->projector.ipImageSize.x, 
						GL_LUMINANCE, GL_FLOAT, &(projector.fImage[0][0])  );		// replace these constants later on

		else
			glDrawPixels(	this->projector.ipImageSize.y, 
							this->projector.ipImageSize.x, 
							GL_RGBA, GL_FLOAT, &(projector.fRGBImage[0][0])  );
		
		glClear(GL_DEPTH_BUFFER_BIT) ;
		this->ShowColorBar() ; 
   }  
}

void VolumeRenderer::SetProjectionResolution(bool bFastPass, bool bIntensity){
	if( bFastPass){			// fast pass projection has been requested
		this->projector.bFastPassMethod = true ;			// subsampling
		this->volumeImg.fRayCastStep = 4 ;					// less dense sampling
		this->bIntensity = true ;							// definitely intensity image
	}
	else{
		this->projector.bFastPassMethod = false ;			// full sampling
		this->volumeImg.fRayCastStep = 1 ;					// dense sampling
		this->bIntensity = bIntensity ;						// can be intensity or color image
	}
	this->bRayCastProjChanged = true ;						// projection changed!
}

void VolumeRenderer::BuildColorBar(){

	for( int i = 0 ; i < ColorBarHeight ; i++){
		float intensity = i*65535.0/ColorBarHeight ;
		float opacity = this->projector.opacityFunction.getValue(intensity) ;  		
		__m128 pixelRGB = this->projector.colorTranFun.getColorValue( intensity) ;
		for( int j = 0 ;  j < ColorBarWidth ; j++){
			fRGBColorBar[i][j][0] = ((float*)&pixelRGB)[0] ; 
			fRGBColorBar[i][j][1] = ((float*)&pixelRGB)[1] ; 
			fRGBColorBar[i][j][2] = ((float*)&pixelRGB)[2] ; 
			fRGBColorBar[i][j][3] = 1 ; 
		}
	}
}

void VolumeRenderer::ShowColorBar(){
//	if ( this->bColorBar ) {
		glRasterPos3f(	-1.0 * projectionWidth/2   * this->projector.fDetectorPixelSize ,			// this is basically the position of the near clipping
						-1.0 * projectionHeight/2  * this->projector.fDetectorPixelSize ,				// plane in object coordinate system	
						this->projector.vSourceLocation.VectorZ - this->projector.fSourceToDetector);
		glDrawPixels(ColorBarWidth, ColorBarHeight, GL_RGBA, GL_FLOAT, fRGBColorBar );
//	}
}