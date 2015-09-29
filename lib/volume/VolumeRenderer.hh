/* written by Abedul Haque, December 2009		*/
/* for CS3610 class project						*/
/* VolumeRenderer class file header				*/

#pragma once

#include "../volume/Volume.hh"
#include "../projector/Projector.hh"
#include "../file_ops/Parameters.hh"
#include "../open_gl.hh"
#include <iostream>

#define ColorBarWidth	20
#define ColorBarHeight	100

class VolumeRenderer{
public:	

	Volume volumeImg;							// Encapsulates the 3D volume
	Projector projector;						// projector of the ray-casting
	Parameters params ;							// parameter file for ray-casting and volume specs
	bool bRayCastProjChanged;					// Do we need to update the projection image
												// if this is set true, we must call CreateVolImage before calling RenderVolImage
	bool bVisible;								// Is volume rendering enabled
	bool bIntensity;							// The defult image is intensity image	
	bool bColorBar;

	float fRGBColorBar[ColorBarHeight][ColorBarWidth][4] ;
	//methods
	VolumeRenderer(){
		bRayCastProjChanged = true ;			// so that there is a projection at first														
		bVisible = false ;						// visible
		bIntensity = false;						// color image
		bColorBar = false ;
	}
	~VolumeRenderer(){}
	
private:
	void InitializeProjector(int type) ;

public:
	void Initialize(std::string szp, int width, int type) ;
	void CreateVolImage() ;
	void RenderVolImage() ;
	void RotateObject(float fAngleX, float fAngleY, float fAngleZ) ;
	void ChangeProjectorLocation(float fDelta) ;
	void SetProjectionResolution( bool bFastPass, bool bIntensity) ;
	void BuildColorBar() ;
	void ShowColorBar() ;

	int projectionWidth ;
	int projectionHeight ;

};