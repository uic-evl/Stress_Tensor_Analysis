#pragma once

#include "../open_gl.hh"
#include <cmath>
#include <iostream>
#include <vector>

#include "../vector/Point.hh"
#include "../vector/Vector.hh"

#include "../draw/Shaders.h"

#define MAX_NUMBER_STREAMLINE	50	

#ifndef min
#define min(a,b) (a<b?a:b)
#endif

#ifndef max
#define max(a,b) (a>b?a:b)
#endif

#define SELECTED	1 
#define UNSELECTED	0

class Streamline
{
	public:
		Vector*** vectorField ;				// vector field on which we are going to calculate streamline
		int slices ;						// number of z-axis slices
		int rows ;							// number of rows
		int cols ;							// number of cols
		int iNumberOfStreamlines ;			// total number of streamlines to be generated
		float fSeedPosX ; //= 97.0 ;
		float fSeedPosY ; //= 120.0 ;
		float fSeedPosZ ; //= 97.0 ;
		
		double boxLenX;
		double boxLenY;
		double boxLenZ;
		
		GLuint StreamlineDisplayListIndex ;	// each streamline will be compiled in a display list. The display list index
											// starting point is saved in streamlineDisplayListIndex. So, streamlineDisplayListIndex is
											// is the index for first streamline, streamlineDisplayListIndex+1 is the index for the second
											// streamline and so on.
		GLuint TempSelectedStreamline ;		// this is display list index of he streamline over which left-button 
											// down mouse cursor is hovering. Note, only one streamline can be selected in this way
		GLuint SelectedStreamlines[MAX_NUMBER_STREAMLINE] ;		//permanently selected streamlines. 0 if the streamline is not selected
																// 1 if it is selected.
		bool bReadyToProcess ;				// whether all data is OK and we can start strealine generation
		bool bVisible;						// whether the streamlines are visible or not. This must set to true if you
											// want to render the streamlines
		float streamlineColor[5][3] ;
		int streamlineColorIndex[MAX_NUMBER_STREAMLINE] ;
		int freeColorForSelected[3] ;

		std::vector<GLfloat> points;
		GLuint vbo[MAX_NUMBER_STREAMLINE];
		int sizes[MAX_NUMBER_STREAMLINE];
		
		GLuint shader; // shader for drawing streamlines

		bool bSeedpointSelected ;
		bool bSeedpointMoved ;

		std::vector<Vector*> stream_pts; // integrated points

	Streamline(void) ;
    ~Streamline(void){}

	void SetVectorField( Vector*** vf, int slices, int rows, int cols){
		this->vectorField = vf ;
		this->slices = slices ;
		this->rows = rows ;
		this->cols = cols ;
	}
	void SetSpacing(float dx, float dy, float dz) ;
	void biLinearInterpolate(Vector &v1, Vector &v2, Vector &v3, Vector &v4, float fDelta1, float fDelta2) ;
	void Configure(fPoint SeedPoint, int n) ;
	void DeleteStreamlines() ;
	void CreateStreamlines() ;
	void RenderStreamlines(GLfloat *eigenfloats) ;
	void SelectModeStreamlinesDraw() ;			// draw the stream lines in GL_SELECT mode
	void SelectModeSeedPointDraw() ;			// draw the seedpoint in GL_SELECT mode
	void UpdateStreamlineStatus(int i) ;
	void TemporarySelectedStreamline(int i) ;
	void DrawSeedPoint() ;
	void UpdateSeedPoint() ;
	void MoveSeedPoint(fPoint& fp) ;
	void GetSeedPoint(fPoint& fp) ;
	
	Vector *getStreamlinePoint(int index){
		return stream_pts[index] ;
	}
	int getSize(){
		return stream_pts.size();
	}
	void clear_points(){
		
		stream_pts.clear();
		points.clear();
	}
	
	
private:
	bool EulerStep(Vector CurrPos, Vector& NextPos, float fStepSize) ;
	bool DerivativeAtPosition(Vector vCurrPos, Vector& vDerivative) ;
	bool MidPointMethod(Vector CurrPos, Vector& NextPos, float fStepSize) ; 
	bool RK4Method(Vector CurrPos, Vector& NextPos, float fStepSize) ; 
};
