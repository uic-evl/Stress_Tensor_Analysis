/* written by Abedul Haque, December 2009		*/
/* for CS3610 class project						*/

#pragma once

#include "../streamlines/Streamline.hh"
#include "../open_gl.hh"
#include "../vector/Vector.hh"
#include "../vector/point.hh"

#define MY_PI 3.14

class VectorFieldRenderer{

public:
	Vector ***VectorField;				// storage for the vector field
	Vector ***EigenValues;				// storage for the eigen values
	
	Streamline *streamline;
	
	int iSlices ;						// number of planes along z axis or least changed index
	int iRows ;							// number of rows or second lest changed index
	int iCols ;							// third index or number of cols

	// these 3 numbers indicates which slice we are showing or rendering
	int iVectorFieldSliceNumberX ;		// slice number along x axis. 
	int iVectorFieldSliceNumberY ;		// along y
	int iVectorFieldSliceNumberZ ;		// along z . Initially we show this slice

	bool bSliceXVisible ;
	bool bSliceYVisible ;
	bool bSliceZVisible ;

	bool bVisible ;							// Global switch for visibility of the vector field. 
											// To see a slice through the vector field, this one has to be set to true
											// and appropriate slice number has to be set.
	
	bool rendered;
	
	GLuint index;							// display list base
	GLuint *lists;							// list indices
	
	GLfloat zoomScale;
	bool scale;
	GLuint scaleIndex;
	GLuint *scaleList;
	
	int skip;
	
	double boxLenX;
	double boxLenY;
	double boxLenZ;
	
	// methods
	VectorFieldRenderer(){
		iVectorFieldSliceNumberX=0;		
		iVectorFieldSliceNumberY=0;		
		iVectorFieldSliceNumberZ=0;		

		bSliceXVisible = false ;
		bSliceYVisible = false ;
		bSliceZVisible = false  ;
		bVisible = false	;	

		VectorField = NULL;
	} 
	void SetSpacing(double dx, double dy, double dz);
	~VectorFieldRenderer(){}
	void Configure(Vector*** vf, iPoint volSize, Vector vSpacing ) ;
	void Configure( Vector*** vf, Vector ***ev, Streamline *s, int slices, int rows, int cols);
	//void setSpacing(float x, float y, float z) ;
	void RenderSlice() ;
	void RenderStreamlines(float mag);
	void IncreaseSliceNumber();
	void DecreaseSliceNumber() ;
	void ToggleSliceVisibilityX() ;
	void ToggleSliceVisibilityY() ;
	void ToggleSliceVisibilityZ() ;
	void drawQuadric(bool shift);
	void drawZoomedQuadric();
	void setZoom(GLfloat zoom);
	void setZoom(GLfloat zoom, float mag, bool redo);
	void cleanup();
	void cleanupzoom();
	
private:
	void DrawAlignedGlyph(Vector v, Vector e) ;
	void DrawSuperQuadric(Vector V);
	void DrawGlyph(void) ;
	void storeQuadric(GLfloat verts[20][20][4], int rs, int vs);
};