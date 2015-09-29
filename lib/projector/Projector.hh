/* this file has been taken from MTW project		*/
/* I modified the file as needed. There have been	*/
/* a lot of changes. I didn't change the comments.	*/
/* So, some comments may be misleading.				*/		


#pragma once

#include "../vector/Vector.hh"
#include "../transform/Transform.hh"
#include "../volume/Volume.hh"
#include "../transform/TransferFunction.hh"
#include "../volume/ColorTransferFunction.hh"


#define AccumulateFour(x) ((double)((float*)&x)[0] + (double)((float*)&x)[1] + (double)((float*)&x)[2] + (double)((float*)&x)[3])
#define MaxFour(x) max(max(((float*)&x)[0],((float*)&x)[1]),max(((float*)&x)[2],((float*)&x)[3]))


typedef struct FPinterp
{
    float   fpWeight[2][2];             // this is the weight given to each true value in the interpolation
}   FPinterp;

typedef struct ResamplePoint
{
    u_short i;                          // these are the coordinates of a point on the edge of a fastpass DRR that need to be resampled
    u_short j;                          //
}   ResamplePoint;

class Projector
{
// public member variables
public:
    char    szProjectorName[MAX_PATH];  // This is just a reference name for the projector object -- used in debugging
    Vector  vSourceLocation;            // the location of the x-ray point source in world coordinates (in mm)
    // The central ray is defined as a vector originating at the x-ray source, perpendicular to the detector screen.
    Vector  vCentralRayOrientation;     // the rotation angles (in degrees) of the central ray vs. world coordinates
    Vector  vCR_Image_Intersection;     // the X-Y screen coordinates where the central ray lands (in screen pixel units)
    Vector  vProjectorScale;            // scale of the projection -- this is the pixel size of the detector screen (in mm).
    float   fSourceToDetector;          // the distance (in mm) from the source to the detector screen along the central ray.
    iPoint  ipImageSize;                // size of detector screen in pixels (for allocating the array)
                                        // note that image coord -Y aligns with Projector coordinate +Z and image +X aligns with Projector +X
    float   fDetectorPixelSize;         // size of square detector pixels (in mm)
    Vector  vPrincipalPosition;         // the location (in screen pixel units) of the central spot -- always (0,S2D,0).
    float** fImage;                     // pointer to standard AllocArray2D array.
	//float fRGBImage[512][512][4] ;				// If you want to do a color projection
	float** fRGBImage;
									// You will need to use ColorTransferFunction for that
	//float** fRefreshBuffer	;			// fImage is cleared with fRefreshBuffer image, Added by Abed
	TransferFunction opacityFunction ;	// Opacity transfer function of the projector
	ColorTransferFunction colorTranFun ;
    Vector  vImageMinScrCoord;          // Bounding rectangle for projection of CT volume on screen
    Vector  vImageMaxScrCoord;          // "
    Transform   tWorld2Projector;       // matrix to transform world coordinates to Projector coordinates
    Transform   tProjector2World;       // matrix to transform Projector coordinates to world coordinates
    int     iProjectionCount;           // number of images projected with this object
    int     iLastProjectionTime;        // number of milliseconds to render the last view
    int     iTotalProjectionTime;       // total time spend in RayCast() function (in milliseconds)
    double  dCorrelation;               // correlation value calculated by this->Correlation();
    float fMaxImageValue;               // Maximum pixel value in the image
	float fMaxRedValue ;				// for color projection
	float fMaxGreenValue ;
	float fMaxBlueValue ;
    float fMeanImageValue;              // Mean pixel value
    float fStdDevImageValue;            // Pixel value standard deviation
    
	enum eSSE
    {
        eSSE_unknown,
        eSSE_SSE2,
        eSSE_nonSSE2
    };
    eSSE SSE_Mode;                      // identifies the mode of processing for SIMD code

    // unnecessary variables
    int     iBeadCount; double  dCalibrationError; double  dRmsErrorInMM;     
    Vector  vLastPositionRendered;   Vector  vLastOrientationRendered;
	///////////////////////

    // fast-pass DRR rendering controls (added on 11-24-2007)
    bool    bFastPassMethod;            // if 'true' RayCast() will plot every Nth row and column, then interpolate to get full image
    int     iFastPassSamplingInterval;  // this is the N in the above line; it comes from the parameters file.
    FPinterp*** fpiInterpolation;       // this table controls the interpolation process for the Nth sampling
    ResamplePoint* rppResample;         // this is a list of the points at the edge of the DRR that we need to resample
    int     nResample;                  // this is a count of the above list
    int     iMaxResamplesPossible;      // this is the limit of the resampling process
public:
    Projector(void);
	~Projector(void);
	
    // Configure the x-ray projector in world space
    void Configure(Vector& XRaySourceLocation, 
                   Vector& CentralRayOrientation, 
                   float   SourceToDetectorDistance, 
                   iPoint& DetectorDimensions, 
                   float   DetectorPixelSize,
                   Vector& CentralRaySpotOnScreen);
    // Set the projector screen buffer to all zeros
    void Clear(void);
    // Project the ctvObj onto the fImage buffer of this Projector; image is added to the present contents of fImage
    void RayCast(Volume& ctvObj);
	// RGB image version of the raycast. The RayCast method generates intensity image. This version generates RGB image
	void Normalize() ;
	// But we need to initialize ColorTransferFunction for that
	void ColorRayCast(Volume& ctvObj);
	void NormalizeRGB() ;
    // Return the average number of milliseconds used per rendering operation
    int AverageRenderTime(void);
    // Calculates the world-to-projector and projector-to-world transforms from the position and orientation members
    void CalculateTransformations(void);
    // Calculate the image coordinates of vWorldPoint's shadow.
    // Note that this routine does not actually plot the point in the fImage array.
    //Vector& ProjectPoint(Vector& vWorldPoint);
    // Prepare the interpolation tables for FastPass processing at the specified sampling interval.
    void CalculateFastPassInterpolationTables(void);
    // Clip the histogram of pixel values down to the specified percentile.
    //void ClipToPercentile(float fClippingPercent);
	// saves the current projection image into the refresh buffer
	//void SaveProjectionToRefreshBuffer(void) ;
	void printProjector() ;
	//void LoadProjector( char *fileName ) ;
};
