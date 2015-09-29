/* this file has been taken from MTW project		*/
/* I modified the file as needed. There have been	*/
/* a lot of changes. I didn't change the comments.	*/
/* So, some comments may be misleading.				*/		

// Projector.cpp

// This class implements an x-ray source and a detector screen.
// It provides a conical projection through a CT object onto the screen.
// The Projector must be positioned and oriented in world coordinates.
// The user must provide a Volume object for projection. (The object contains its position and orientation.)
// The user is responsible for post-processing the image array after projection.

#define _USE_MATH_DEFINES

#include "math.h"
#include "assert.h"
#include "Projector.hh"
#include "../std/AllocArray.hh"
#include "emmintrin.h"      // header file for SSE intrinsics
#include "../std/FatalError.hh"

Projector::Projector(void)
{
    // clear out the object initially
    memset(this, 0, sizeof(Projector));
	
	//If the above block is uncommented, this line should be commented
	this->SSE_Mode = this->eSSE_SSE2;
}

Projector::~Projector(void)
{
	if(fImage != NULL)
		free(fImage) ;
	if(fRGBImage != NULL)
		free(fRGBImage) ;
}

// Configure the x-ray projector in world space
void Projector::Configure(Vector&   XRaySourceLocation, 
                          Vector&   CentralRayOrientation, 
                          float     SourceToDetectorDistance, 
                          iPoint&   DetectorDimensions, 
                          float     DetectorPixelSize,
                          Vector&   CentralRaySpotOnScreen)
{
    // allocate the detector screen buffer
    this->ipImageSize               = DetectorDimensions;
    this->fImage                    = (float**)AllocArray2D(ipImageSize.y,
                                                            ipImageSize.x,
                                                            sizeof(float));
 	this->fRGBImage					= (float **)AllocArray2D(ipImageSize.y, ipImageSize.x*4, sizeof(float)) ;

    this->vCR_Image_Intersection    = CentralRaySpotOnScreen;
    // capture the remaining configuration parameters
    this->vSourceLocation           = XRaySourceLocation  ;
	//this->vSourceLocation.ScalarMult(-1) ;									// -1 is to get the negative distance. So that a translation takes you to origin
    this->vCentralRayOrientation    = CentralRayOrientation;
    this->fSourceToDetector         = SourceToDetectorDistance;
    this->vProjectorScale           = new Vector(DetectorPixelSize,
                                             DetectorPixelSize,
                                             DetectorPixelSize);
    this->fDetectorPixelSize        = DetectorPixelSize;
    this->vPrincipalPosition        = new Vector(0.0f,
                                             0.0f,
                                             -SourceToDetectorDistance / DetectorPixelSize);		// location of detector in projector coordinate
																									// since projector projects towards negative z axis

    // calculate the world-to-Projector transformation
    CalculateTransformations();
    // clear out the performance counters
    this->iProjectionCount = 0;
    this->iTotalProjectionTime = 0;
}

// Set the projector screen buffer to all zeros
void Projector::Clear(void)
{
    memset(&fImage[0][0], 
           0, 
           ipImageSize.x * ipImageSize.y * sizeof(fImage[0][0]) );

    memset(&fRGBImage[0][0], 
           0, 
           ipImageSize.x * ipImageSize.y * 4 * sizeof(fRGBImage[0][0]) );

    this->fMaxImageValue = 0.0;
	this->fMaxRedValue = 0.0 ;
	this->fMaxGreenValue = 0.0 ;
	this->fMaxBlueValue  = 0.0 ;
}

// Project the ctvObj onto the fImage buffer of this Projector; image is added to the present contents of fImage
void Projector::RayCast(Volume& ctvObj)
{
    // set up the FastPass interpolation tables on the first time through
    if (this->fpiInterpolation == NULL) CalculateFastPassInterpolationTables();

    // set up the pixel step values for fast-pass method vs. normal method
    int jStep                   = 1,    // these are the step values for the normal method.
        iStep                   = 1,    //
        iSimdStep               = 4;    //
    if (this->bFastPassMethod)
    {	
        jStep       = iFastPassSamplingInterval;        // these are the step values for the 'fast-pass' method
        iStep       = iFastPassSamplingInterval;        
        iSimdStep   = iFastPassSamplingInterval * 4;    

        // we must allocate the resample list, if it has not already been done
        if (this->rppResample == NULL)
        {
            this->iMaxResamplesPossible = 2 * (this->ipImageSize.x + this->ipImageSize.y) * (this->iFastPassSamplingInterval + 1);
            this->rppResample = (ResamplePoint*) calloc(iMaxResamplesPossible, sizeof(ResamplePoint));
        }   // if (this->rppResample == NULL)

    }   // if (this->bFastPassMethod)

    // set up accumulators for the image
    // calculate the two required transformation matrices
    static Transform tObj2Proj,
                     tProj2Obj;
    ctvObj.tCtVolToWorld.LeftMult(this->tWorld2Projector, tObj2Proj);   // make the transform to go from CT space to projector space
    tProj2Obj.Invert(tObj2Proj);                                        // invert it to go from projector space to CT space.
    // map the vertices of the object to projector space and find the limits;
    // these represent the projector frustum coordinates in Projector space
    static Vector vImageMin;
    static Vector vImageMax;
    vImageMin = new Vector( 1e20f, 1e20f, 1e20f);
    vImageMax = new Vector(-1e20f,-1e20f,-1e20f);
    static Vector vObjVertexInCtSpace   (0,0,0),
                  vObjVertexInProjSpace (0,0,0),
                  vObjVertexOnScreen    (0,0,0);

#define FindFrustum(fVertexX,fVertexY,fVertexZ)                                                                      \
    {                                                                                                                \
        vObjVertexInCtSpace.VectorX = fVertexX;                                                                      \
        vObjVertexInCtSpace.VectorY = fVertexY;                                                                      \
        vObjVertexInCtSpace.VectorZ = fVertexZ;                                                                      \
        tObj2Proj.RightMult(vObjVertexInCtSpace,vObjVertexInProjSpace);												 \
		vObjVertexOnScreen.VectorX = vObjVertexInProjSpace.X() * vPrincipalPosition.Z() / vObjVertexInProjSpace.Z(); \
        vObjVertexOnScreen.VectorZ = vPrincipalPosition.Z();														 \
        vObjVertexOnScreen.VectorY = vObjVertexInProjSpace.Y() * vPrincipalPosition.Z() / vObjVertexInProjSpace.Z(); \
        vImageMin.Min(vObjVertexOnScreen);                                                                           \
        vImageMax.Max(vObjVertexOnScreen);                                                                           \
    }

    FindFrustum(ctvObj.vMinNonzeroPlanes.X(), ctvObj.vMinNonzeroPlanes.Y(), ctvObj.vMinNonzeroPlanes.Z());
    FindFrustum(ctvObj.vMaxNonzeroPlanes.X(), ctvObj.vMinNonzeroPlanes.Y(), ctvObj.vMinNonzeroPlanes.Z());
    FindFrustum(ctvObj.vMinNonzeroPlanes.X(), ctvObj.vMaxNonzeroPlanes.Y(), ctvObj.vMinNonzeroPlanes.Z());
    FindFrustum(ctvObj.vMaxNonzeroPlanes.X(), ctvObj.vMaxNonzeroPlanes.Y(), ctvObj.vMinNonzeroPlanes.Z());
    FindFrustum(ctvObj.vMinNonzeroPlanes.X(), ctvObj.vMinNonzeroPlanes.Y(), ctvObj.vMaxNonzeroPlanes.Z());
    FindFrustum(ctvObj.vMaxNonzeroPlanes.X(), ctvObj.vMinNonzeroPlanes.Y(), ctvObj.vMaxNonzeroPlanes.Z());
    FindFrustum(ctvObj.vMinNonzeroPlanes.X(), ctvObj.vMaxNonzeroPlanes.Y(), ctvObj.vMaxNonzeroPlanes.Z());
    FindFrustum(ctvObj.vMaxNonzeroPlanes.X(), ctvObj.vMaxNonzeroPlanes.Y(), ctvObj.vMaxNonzeroPlanes.Z());

#undef FindFrustum
    // all of the above calculations were done in Projector space.  Now we need to calculate the minimum and maximum
    // in screen coordinates and convert to integer values:

    vImageMinScrCoord = new Vector(vImageMin.X() + vCR_Image_Intersection.X(),  // Image X is Projector X, offset by central spot X
                               vImageMin.Y() + vCR_Image_Intersection.Y(),  // Image Y is Projector Z, offset by central spot Y
                               0.0f);                                       // Image Z is irrelevant
    vImageMinScrCoord.Floor();
    vImageMaxScrCoord = new Vector(vImageMax.X() + vCR_Image_Intersection.X(),
                               vImageMax.Y() + vCR_Image_Intersection.Y(),  
                               0.0f);                                       
    vImageMaxScrCoord.Ceil();
    // now we will clip the screen coordinates to the size specified by ipImageSize
    vImageMinScrCoord.VectorX = max(vImageMinScrCoord.VectorX, 0.0f);
    vImageMinScrCoord.VectorY = max(vImageMinScrCoord.VectorY, 0.0f);
    vImageMaxScrCoord.VectorX = min(vImageMaxScrCoord.VectorX, (float)(ipImageSize.x - 1));
    vImageMaxScrCoord.VectorY = min(vImageMaxScrCoord.VectorY, (float)(ipImageSize.y - 1));

    // declare variables to define the frustum limits
    int jStart, jEnd, iStart, iEnd;
    // use one of two raycasting processes, depending on the hardware in THIS machine
   
        // set up accumulators for statistics
        static __m128 m128_Sum, m128_SumSqr, m128_Max;
        m128_Sum        = _mm_setzero_ps();
        m128_SumSqr     = _mm_setzero_ps();
        m128_Max        = _mm_setzero_ps();
        static int iTupleCnt;
        iTupleCnt       = 0;
        // locate the x-ray source in CT object space
        static Vector  vNull(0,0,0),   // source in projector coordinates; always (0,0,0)
                       vSource;        // source in CT coordinates
        tProj2Obj.RightMult(vNull,vSource);
        // now we can trace each ray through the frustum and add the footprint on the screen
        jStart = (int)vImageMinScrCoord.VectorY,
        jEnd   = (int)vImageMaxScrCoord.VectorY,
        iStart = (int)vImageMinScrCoord.VectorX,
        iEnd   = (int)vImageMaxScrCoord.VectorX;
        for (int jScreen=jStart; jScreen<=jEnd; jScreen+=jStep)         // note that we are stepping wider if bFastPassMethod is 'true'.
            for (int iScreen=iStart; iScreen<=iEnd; iScreen+=iSimdStep) // process four pixels at a whack in the SIMD variables
                                                                        // but skip more X-spaces if bFastPassMethod is 'true'.
        {
            // this is the accumulator that will gather the values of the four pixels
            static __m128 m128_Accum;
			static __m128 m128_Intensity, m128_Opacity ;

            m128_Accum = _mm_set_ps1(0.0f);
            // first get the locations in Projector space of four X-adjacent pixels
			static Vector4 vPixelInCtObjCoords, vPixelInProjectorCoords;
			
			vPixelInProjectorCoords = new	   Vector4((0.0f),                                        // this coordinate is filled in below
											  (float)(jScreen - vCR_Image_Intersection.Y()), //Projector Y is image Y  	
											  (float)(vPrincipalPosition.Z()));                // source-to-detector in Projector units
                                             
            for (int i=0; i<4; i++) // load four x-coordinates
            {
                // set the x coordinate of the four pixels
                int ii = iScreen + i * iStep;   // note that we are stepping by two if bFastPassMethod is 'true'
                ((float*)&vPixelInProjectorCoords.VectorX)[i] = (float)(ii - vCR_Image_Intersection.X());// Projector X is minus Image X, 
                                                                                                        // offset by the central spot
				//vPixelInProjectorCoords.VectorX)[i] = (float)(ii - vCR_Image_Intersection.X());// Projector X is minus Image X, 
            }
            // transform the 4 pixels to CT object coordinates
            tProj2Obj.RightMult(vPixelInProjectorCoords,vPixelInCtObjCoords);
            // compute the sampling delta Vector4
            vPixelInCtObjCoords.Subtract(vSource);  // this is the net ray Vector4
            static __m128 m128_RayLength, m128_Normalizer;
            // calculate the magnitudes of the four rays
			// _mm_add_ps  only adds two m128 variables. So, we need to use it twice to add three variables, x, y and z
            m128_RayLength = _mm_sqrt_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(vPixelInCtObjCoords.VectorX,vPixelInCtObjCoords.VectorX),
                                                               _mm_mul_ps(vPixelInCtObjCoords.VectorY,vPixelInCtObjCoords.VectorY)),
                                                               _mm_mul_ps(vPixelInCtObjCoords.VectorZ,vPixelInCtObjCoords.VectorZ) ));
            // calculate a factor to normalize the four rays to unit length (ctvObj.fRayCastStep)
            m128_Normalizer = _mm_div_ps(_mm_set_ps1(ctvObj.fRayCastStep),m128_RayLength);
            // now normalize the four rays to unit length (ctvObj.fRayCastStep)
            vPixelInCtObjCoords.ScalarMult(&m128_Normalizer); // vPixelInCtObjCoords has been converted into the step delta Vector4
#define vRayCastStepDelta   vPixelInCtObjCoords      // we will now use a different name for the same variable
            // Add a 'Kahan zero' to the step delta to avoid a possible divide-by-zero
            vRayCastStepDelta.VectorX = _mm_add_ps(vRayCastStepDelta.VectorX,_mm_set_ps1(1e-20f));
            vRayCastStepDelta.VectorY = _mm_add_ps(vRayCastStepDelta.VectorY,_mm_set_ps1(1e-20f));
            vRayCastStepDelta.VectorZ = _mm_add_ps(vRayCastStepDelta.VectorZ,_mm_set_ps1(1e-20f));
            // find the intersection points of the rays with the bounding planes
            // all of these values are in units of "number of steps"
            static __m128 m128_nx1,m128_nx2,m128_nxL,m128_nxU;
            static __m128 m128_ny1,m128_ny2,m128_nyL,m128_nyU;
            static __m128 m128_nz1,m128_nz2,m128_nzL,m128_nzU;
            {
                m128_nx1 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMinNonzeroPlanes.VectorX),  // x intersection 1
                                                 _mm_set_ps1(vSource.VectorX)                  ),
                                        vRayCastStepDelta.VectorX);
                m128_nx2 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMaxNonzeroPlanes.VectorX),  // x intersection 2
                                                 _mm_set_ps1(vSource.VectorX)                  ),
                                        vRayCastStepDelta.VectorX);
                m128_nxL = _mm_min_ps(m128_nx1,m128_nx2);                                        // first x intersection
                m128_nxU = _mm_max_ps(m128_nx1,m128_nx2);                                        // last x intersection
                m128_ny1 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMinNonzeroPlanes.VectorY),  // y intersection 1
                                                 _mm_set_ps1(vSource.VectorY)                  ),
                                        vRayCastStepDelta.VectorY);
                m128_ny2 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMaxNonzeroPlanes.VectorY),  // y intersection 2
                                                 _mm_set_ps1(vSource.VectorY)                  ),
                                        vRayCastStepDelta.VectorY);
                m128_nyL = _mm_min_ps(m128_ny1,m128_ny2);                                        // first y intersection
                m128_nyU = _mm_max_ps(m128_ny1,m128_ny2);                                        // last y intersection
                m128_nz1 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMinNonzeroPlanes.VectorZ),  // z intersection 1
                                                 _mm_set_ps1(vSource.VectorZ)                  ),
                                        vRayCastStepDelta.VectorZ);
                m128_nz2 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMaxNonzeroPlanes.VectorZ),  // z intersection 2
                                                 _mm_set_ps1(vSource.VectorZ)                  ),
                                        vRayCastStepDelta.VectorZ);
                m128_nzL = _mm_min_ps(m128_nz1,m128_nz2);                                        // first z intersection
                m128_nzU = _mm_max_ps(m128_nz1,m128_nz2);                                        // last z intersection
            }
            // find the start and end points on the four rays, considering the bounding planes
            static __m128 m128_nStart,m128_nEnd;
            static __m128 m128_Count;
            m128_nStart  = _mm_max_ps(m128_nxL,_mm_max_ps(m128_nyL,m128_nzL));    // point of piercing the last of the closest planes (enter)
            m128_nEnd    = _mm_min_ps(m128_nxU,_mm_min_ps(m128_nyU,m128_nzU));    // point of piercing the first of the furthest planes (exit)
            m128_Count   = _mm_sub_ps(m128_nEnd,m128_nStart);                     // number of samples required on each ray
            // if the ray does not pierce the bounded CT object, the count will be <= 0
            int iSampleCount = 0;
            for (int i=0; i<4; i++)
            {
                // calculate the maximum number of samples we must take in order to get the longest ray
                iSampleCount = max(iSampleCount, (int)((float*)&m128_Count)[i]);  // drop fractional part
            }
            // calculate the first sample point for each ray
            static Vector4 vSamplePoint;
            vSamplePoint = vRayCastStepDelta;
            if (iSampleCount <= 0) goto NoPiercing; // we can skip these four rays completely if none pierce the CT object!
            vSamplePoint.ScalarMult(&m128_nStart);
            vSamplePoint.Add(vSource);
            // now we can run the rays
            for (int iSample=0; iSample<=iSampleCount; iSample++)
            {
                // round the indices to nearest neighbor
                static __m128i m128i_i,m128i_j,m128i_k;
                m128i_i = _mm_cvtps_epi32(vSamplePoint.VectorX);
                m128i_j = _mm_cvtps_epi32(vSamplePoint.VectorY);
                m128i_k = _mm_cvtps_epi32(vSamplePoint.VectorZ);
                // accumulate an __m128 of voxels by picking the four values out of the CT object
                for (int i=0; i<4; i++)
                {
                    if (((float*)&m128_Count)[i] >= iSample)
                    {
                        assert(((int*)&m128i_k)[i] >= ctvObj.vMinNonzeroPlanes.VectorZ && 
                               ((int*)&m128i_k)[i] <= ctvObj.vMaxNonzeroPlanes.VectorZ);
                        assert(((int*)&m128i_j)[i] >= ctvObj.vMinNonzeroPlanes.VectorY && 
                               ((int*)&m128i_j)[i] <= ctvObj.vMaxNonzeroPlanes.VectorY);
                        assert(((int*)&m128i_i)[i] >= ctvObj.vMinNonzeroPlanes.VectorX && 
                               ((int*)&m128i_i)[i] <= ctvObj.vMaxNonzeroPlanes.VectorX);


                        ((float*)&m128_Accum)[i] +=  this->opacityFunction.getValue((float)ctvObj.usCtVol[((int*)&m128i_k)[i]]
                                                                        [((int*)&m128i_j)[i]]
                                                                        [((int*)&m128i_i)[i]] );

                    }

					
                }
                // index to the next four sample positions
                vSamplePoint.VectorX = _mm_add_ps(vSamplePoint.VectorX,vRayCastStepDelta.VectorX);
                vSamplePoint.VectorY = _mm_add_ps(vSamplePoint.VectorY,vRayCastStepDelta.VectorY);
                vSamplePoint.VectorZ = _mm_add_ps(vSamplePoint.VectorZ,vRayCastStepDelta.VectorZ);
            }

            // put away the rays differently, depending on bFastPassMethod switch
            if (bFastPassMethod)
            {
                fImage[jScreen][iScreen                                ] = ((float*)&m128_Accum)[0];
                fImage[jScreen][iScreen +     iFastPassSamplingInterval] = ((float*)&m128_Accum)[1];
                fImage[jScreen][iScreen + 2 * iFastPassSamplingInterval] = ((float*)&m128_Accum)[2];
                fImage[jScreen][iScreen + 3 * iFastPassSamplingInterval] = ((float*)&m128_Accum)[3];
            }   // if (this->bFastPassMethod)
            else    // if ( ! this->bFastPassMethod) //i.e., normal method
            {
                fImage[jScreen][iScreen    ] = ((float*)&m128_Accum)[0];
                fImage[jScreen][iScreen + 1] = ((float*)&m128_Accum)[1];
                fImage[jScreen][iScreen + 2] = ((float*)&m128_Accum)[2];
                fImage[jScreen][iScreen + 3] = ((float*)&m128_Accum)[3];
            }   // else    // if ( ! this->bFastPassMethod) //i.e., normal method

            // accumulate the stats
            m128_Sum    = _mm_add_ps(m128_Sum,m128_Accum);
            m128_SumSqr = _mm_add_ps(m128_SumSqr,_mm_mul_ps(m128_Accum,m128_Accum));
            m128_Max    = _mm_max_ps(m128_Max,m128_Accum);
            iTupleCnt++;
NoPiercing:
            continue;
        }   // for (int jScreen=jStart; jScreen<=jEnd; jScreen+=jStep) for (int iScreen=iStart; iScreen<=iEnd; iScreen+=iSimdStep)
        // finish the stats
        fMeanImageValue          = (float)(AccumulateFour(m128_Sum)    / (4.0f * iTupleCnt));
        float fMeanSqrImageValue = (float)(AccumulateFour(m128_SumSqr) / (4.0f * iTupleCnt));
        fStdDevImageValue        = sqrtf(fMeanSqrImageValue - sqr(fMeanImageValue));
        fMaxImageValue           = MaxFour(m128_Max);

        // if we are in FastPass mode, this is where we refine the edges of the DRR (AKA "The Stephanie Edge")
        if (this->bFastPassMethod)
        {
            // go back and recheck every point we projected to find the actual edge of the DRR
            this->nResample = 0;
            for (int jScreen=jStart; jScreen<=jEnd; jScreen+=jStep) 
                for (int iScreen=iStart; iScreen<=iEnd; iScreen+=iStep)
            {
                // check to be sure that we are not past the edge; if we are, adjust it
                int jInternalStep = min(jStep, ipImageSize.y - jScreen - 1);
                int iInternalStep = min(iStep, ipImageSize.x - iScreen - 1);
                /*
                if ((jScreen + jStep) >= this->ipImageSize.y || (iScreen + iStep) >= this->ipImageSize.x) continue;
                */

                // now check the four corners of the interpolation region to see how many are there
                int iCornerCount = 0;
                if (fImage[jScreen                ][iScreen                ] > 0) iCornerCount++;
                if (fImage[jScreen                ][iScreen + iInternalStep] > 0) iCornerCount++;
                if (fImage[jScreen + jInternalStep][iScreen                ] > 0) iCornerCount++;
                if (fImage[jScreen + jInternalStep][iScreen + iInternalStep] > 0) iCornerCount++;
                if (iCornerCount == 0)  continue;   // there is no data here whatsoever
                if (iCornerCount == 4)  continue;   // there is enough data to provide a good interpolation; we are in the middle

                // if we have found 1, 2, or 3 corners, this square is on the edge and needs to be finely sampled
                for (int j=jScreen; j<=jScreen+jInternalStep; j++) 
                    for (int i=iScreen; i<=iScreen+iInternalStep; i++) 
                {
                    // first check to see if this pixel is already covered
                    if (fImage[j][i] != 0) continue;    // this pixel is already considered

                    // now put it on the list and mark it so we don't consider it again
                    this->rppResample[this->nResample].i = (u_short) i;
                    this->rppResample[this->nResample].j = (u_short) j;
                    this->nResample++;
                    if (this->nResample >= iMaxResamplesPossible) break;    // don't overrun our buffer
                    fImage[j][i] = -1.0e-30f;   // mark the pixel with a tiny negative value so it won't affect any calculations

                }   // for (int j=jScreen; j<=jScreen+jInternalStep; j++) for (int i=iScreen; i<=iScreen+iInternalStep; i++)

            }   // for (int jScreen=jStart; jScreen<=jEnd; jScreen+=jStep) for (int iScreen=iStart; iScreen<=iEnd; iScreen+=iStep)

            // we now have a list (rppResample[]) of all the edge point that need to be resampled
            // we will process them in groups of four pixels so we can use SIMD arithmetic
            for (int iResampleIndex = 0; iResampleIndex < this->nResample; iResampleIndex += 4)
            {
                // this is the accumulator that will gather the values of the four pixels
                static __m128 m128_Accum;
                m128_Accum = _mm_set_ps1(0.0f);

                // first get the locations in Projector space of four X-adjacent pixels
                static Vector4 vPixelInProjectorCoords,
                            vPixelInCtObjCoords;

				vPixelInProjectorCoords = new Vector4((0.0f),                                        // this coordinate is filled in below
											  (float)(0.0f),	
											  (float)(vPrincipalPosition.Z()));                // source-to-detector in Projector units

                for (int i=0; i<4; i++) // load four sets of coordinates
                {
                    int iInnerLoopIndex = iResampleIndex + i;

                    if (iInnerLoopIndex >= this->nResample)
                    {
                        // this is NOT a real pixel of this resample; let's nail it to the image origin so it is out of the way
                        memset(&rppResample[iInnerLoopIndex], 0, sizeof(ResamplePoint));
                    }

                    int ii = rppResample[iInnerLoopIndex].i;                                                    // Image X coordinate of pixel to be sampled
                    ((float*)&vPixelInProjectorCoords.VectorX)[i] = (float)(ii - vCR_Image_Intersection.X());    // Projector X is minus Image X, offset by the central spot
                        
                    int jj = rppResample[iInnerLoopIndex].j;                                                    // Image Y coordinate of pixel to be sampled
                    ((float*)&vPixelInProjectorCoords.VectorY)[i] =  (float)(jj - vCR_Image_Intersection.Y());    // Projector Z is Image Y, offset by the central spot

                }

                // transform the 4 pixels to CT object coordinates
                tProj2Obj.RightMult(vPixelInProjectorCoords,vPixelInCtObjCoords);

                // compute the sampling delta Vector4
                vPixelInCtObjCoords.Subtract(vSource);  // this is the net ray Vector4
                static __m128 m128_RayLength, m128_Normalizer;

                // calculate the magnitudes of the four rays
                m128_RayLength = _mm_sqrt_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(vPixelInCtObjCoords.VectorX,vPixelInCtObjCoords.VectorX),
                                                                   _mm_mul_ps(vPixelInCtObjCoords.VectorY,vPixelInCtObjCoords.VectorY)),
                                                                   _mm_mul_ps(vPixelInCtObjCoords.VectorZ,vPixelInCtObjCoords.VectorZ) ));

                // calculate a factor to normalize the four rays to unit length (ctvObj.fRayCastStep)
                m128_Normalizer = _mm_div_ps(_mm_set_ps1(ctvObj.fRayCastStep),m128_RayLength);

                // now normalize the four rays to unit length (ctvObj.fRayCastStep)
                vPixelInCtObjCoords.ScalarMult(&m128_Normalizer); // vPixelInCtObjCoords has been converted into the step delta Vector4

                // Add a 'Kahan zero' to the step delta to avoid a possible divide-by-zero
                vRayCastStepDelta.VectorX = _mm_add_ps(vRayCastStepDelta.VectorX,_mm_set_ps1(1e-20f));
                vRayCastStepDelta.VectorY = _mm_add_ps(vRayCastStepDelta.VectorY,_mm_set_ps1(1e-20f));
                vRayCastStepDelta.VectorZ = _mm_add_ps(vRayCastStepDelta.VectorZ,_mm_set_ps1(1e-20f));

                // find the intersection points of the rays with the bounding planes
                // all of these values are in units of "number of steps"
                static __m128 m128_nx1,m128_nx2,m128_nxL,m128_nxU;
                static __m128 m128_ny1,m128_ny2,m128_nyL,m128_nyU;
                static __m128 m128_nz1,m128_nz2,m128_nzL,m128_nzU;
                {
                    m128_nx1 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMinNonzeroPlanes.VectorX),  // x intersection 1
                                                     _mm_set_ps1(vSource.VectorX)                  ),
                                          vRayCastStepDelta.VectorX);
                    m128_nx2 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMaxNonzeroPlanes.VectorX),  // x intersection 2
                                                     _mm_set_ps1(vSource.VectorX)                  ),
                                          vRayCastStepDelta.VectorX);
                    m128_nxL = _mm_min_ps(m128_nx1,m128_nx2);                                        // first x intersection
                    m128_nxU = _mm_max_ps(m128_nx1,m128_nx2);                                        // last x intersection
                    m128_ny1 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMinNonzeroPlanes.VectorY),  // y intersection 1
                                                     _mm_set_ps1(vSource.VectorY)                  ),
                                          vRayCastStepDelta.VectorY);
                    m128_ny2 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMaxNonzeroPlanes.VectorY),  // y intersection 2
                                                     _mm_set_ps1(vSource.VectorY)                  ),
                                          vRayCastStepDelta.VectorY);
                    m128_nyL = _mm_min_ps(m128_ny1,m128_ny2);                                        // first y intersection
                    m128_nyU = _mm_max_ps(m128_ny1,m128_ny2);                                        // last y intersection
                    m128_nz1 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMinNonzeroPlanes.VectorZ),  // z intersection 1
                                                     _mm_set_ps1(vSource.VectorZ)                  ),
                                          vRayCastStepDelta.VectorZ);
                    m128_nz2 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMaxNonzeroPlanes.VectorZ),  // z intersection 2
                                                     _mm_set_ps1(vSource.VectorZ)                  ),
                                          vRayCastStepDelta.VectorZ);
                    m128_nzL = _mm_min_ps(m128_nz1,m128_nz2);                                        // first z intersection
                    m128_nzU = _mm_max_ps(m128_nz1,m128_nz2);                                        // last z intersection
                }

                // find the start and end points on the four rays, considering the bounding planes
                static __m128 m128_nStart,m128_nEnd;
                static __m128 m128_Count;
                m128_nStart  = _mm_max_ps(m128_nxL,_mm_max_ps(m128_nyL,m128_nzL));    // point of piercing the last of the closest planes (enter)
                m128_nEnd    = _mm_min_ps(m128_nxU,_mm_min_ps(m128_nyU,m128_nzU));    // point of piercing the first of the furthest planes (exit)
                m128_Count   = _mm_sub_ps(m128_nEnd,m128_nStart);                     // number of samples required on each ray

                // if the ray does not pierce the bounded CT object, the count will be <= 0
                int iSampleCount = 0;
                for (int i=0; i<4; i++)
                {
                    // calculate the maximum number of samples we must take in order to get the longest ray
                    iSampleCount = max(iSampleCount, (int)((float*)&m128_Count)[i]);  // drop fractional part
                }

                // calculate the first sample point for each ray
                static Vector4 vSamplePoint;
                vSamplePoint = vRayCastStepDelta;
                if (iSampleCount <= 0) goto NoPiercing2; // we can skip these four rays completely if none pierce the CT object!
                vSamplePoint.ScalarMult(&m128_nStart);
                vSamplePoint.Add(vSource);

                // now we can run the rays
                for (int iSample=0; iSample<=iSampleCount; iSample++)
                {
                    // round the indices to nearest neighbor
                    static __m128i m128i_i,m128i_j,m128i_k;
                    m128i_i = _mm_cvtps_epi32(vSamplePoint.VectorX);
                    m128i_j = _mm_cvtps_epi32(vSamplePoint.VectorY);
                    m128i_k = _mm_cvtps_epi32(vSamplePoint.VectorZ);

                    // accumulate an __m128 of voxels by picking the four values out of the CT object
                    for (int i=0; i<4; i++)
                    {
                        if (((float*)&m128_Count)[i] >= iSample)
                        {
                            assert( ((int*)&m128i_k)[i] >= ctvObj.vMinNonzeroPlanes.VectorZ && 
                                    ((int*)&m128i_k)[i] <= ctvObj.vMaxNonzeroPlanes.VectorZ);
                            assert( ((int*)&m128i_j)[i] >= ctvObj.vMinNonzeroPlanes.VectorY && 
                                    ((int*)&m128i_j)[i] <= ctvObj.vMaxNonzeroPlanes.VectorY);
                            assert( ((int*)&m128i_i)[i] >= ctvObj.vMinNonzeroPlanes.VectorX && 
                                    ((int*)&m128i_i)[i] <= ctvObj.vMaxNonzeroPlanes.VectorX);
                            //m128_Accum)[i] +=  (float)ctvObj.usCtVol[m128i_k)[i]]
                            //                                                [m128i_j)[i]]
                            //                                                [m128i_i)[i]];

							((float*)&m128_Accum)[i] +=  this->opacityFunction.getValue((float)ctvObj.usCtVol[((int*)&m128i_k)[i]]
																			[((int*)&m128i_j)[i]]
	                                                                        [((int*)&m128i_i)[i]] );

                        }
                    }

                    // index to the next four sample positions
                    vSamplePoint.VectorX = _mm_add_ps(vSamplePoint.VectorX,vRayCastStepDelta.VectorX);
                    vSamplePoint.VectorY = _mm_add_ps(vSamplePoint.VectorY,vRayCastStepDelta.VectorY);
                    vSamplePoint.VectorZ = _mm_add_ps(vSamplePoint.VectorZ,vRayCastStepDelta.VectorZ);

                }

                for (int i=0; i<4; i++) // now store the results at the four sets of coordinates
                {
                    fImage[rppResample[iResampleIndex + i].j][rppResample[iResampleIndex + i].i] = ((float*)&m128_Accum)[i];
                }

NoPiercing2:    continue;

            }   // for (int iResampleIndex = 0; iResampleIndex < this->nResample; iResampleIndex += 4)

        }   // if (this->bFastPassMethod)

#ifdef WriteProjectionFile
    {
        char szFileName[MAX_PATH];
        sprintf(szFileName,
                "%sDRRbeforeInterpolation.raw",
                this->szProjectorName);
        FILE* fpImageFile = fopen(szFileName, "wb");
        if (fpImageFile)
        {
            fwrite( &fImage[0][0], 
                    sizeof(float), 
                    ipImageSize.x * ipImageSize.y, 
                    fpImageFile);
            fclose( fpImageFile);
        }
        else
        {
            FATAL_ERROR(szFileName);
        }
    }   // #ifdef WriteProjectionFile
#endif 

    if (this->bFastPassMethod)  // we must now interpolate the partially projected CT scan into a full-density image
    {
        for (int jScreen=jStart; jScreen<=(jEnd-jStep); jScreen+=jStep) 
            for (int iScreen=iStart; iScreen<=(iEnd-iStep); iScreen+=iStep)
        {
            // pull out the 'real' sampled values surrounding the interpolation region
            float fSample00 = fImage[jScreen        ][iScreen        ];
            float fSample01 = fImage[jScreen        ][iScreen + iStep];
            float fSample10 = fImage[jScreen + jStep][iScreen        ];
            float fSample11 = fImage[jScreen + jStep][iScreen + iStep];

            // determine the interpolation mode (see CalculateFastPassInterpolationTables())
            int iInterpMode = 4;
            if (fSample00 == 0 ||
                fSample01 == 0 ||
                fSample10 == 0 ||
                fSample11 == 0 )    continue;   // we must have all four corner posts in order to interpolate

            // now spin through the interpolated pixels
            for (int jInter=0; jInter<=iFastPassSamplingInterval; jInter++) 
                for (int iInter=0; iInter<=iFastPassSamplingInterval; iInter++) 
            {
                // calculate the interpolated pixel indices
                int ii = iScreen + iInter;
                int jj = jScreen + jInter;

                // check to see if this pixel has already been calculated
                if (fImage[jj][ii]) continue;   // no need to recalculate one that's already been done, such as a sampled pixel

                // perform the interpolation
                fImage[jj][ii] =    fSample00 * fpiInterpolation[iInterpMode][jInter][iInter].fpWeight[0][0] +
                                    fSample01 * fpiInterpolation[iInterpMode][jInter][iInter].fpWeight[0][1] +
                                    fSample10 * fpiInterpolation[iInterpMode][jInter][iInter].fpWeight[1][0] +
                                    fSample11 * fpiInterpolation[iInterpMode][jInter][iInter].fpWeight[1][1] ;

            }   // for (int jInter=0; jInter<=iFastPassSamplingInterval; jInterp++) for (int iInter=0; iInter<=iFastPassSamplingInterval; iInterp++) 

        }   // for (int jScreen=jStart; jScreen<=jEnd; jScreen+=jStep) for (int iScreen=(iStart+1); iScreen<=iEnd; iScreen+=iStep)

    }   // if (this->bFastPassMethod)  // we must now interpolate the partially projected CT scan into a full-density image

#ifdef WriteProjectionFile
    {
        char szFileName[MAX_PATH];
        sprintf(szFileName,
                "%sDRRafterInterpolation.raw",
                this->szProjectorName);
        FILE* fpImageFile = fopen(szFileName, "wb");
        if (fpImageFile)
        {
            fwrite( &fImage[0][0], 
                    sizeof(float), 
                    ipImageSize.x * ipImageSize.y, 
                    fpImageFile);
            fclose( fpImageFile);
        }
        else
        {
            FATAL_ERROR(szFileName);
        }
    }   // #ifdef WriteProjectionFile
#endif 

    // mark the end time and accumulate it

	// now we will extend the boundary if there was any more projection already in the frameBuffer

    this->iProjectionCount++;

}

void Projector::Normalize(){
        int jStart = (int)vImageMinScrCoord.VectorY ;
        int jEnd   = (int)vImageMaxScrCoord.VectorY ;
        int iStart = (int)vImageMinScrCoord.VectorX ;
        int iEnd   = (int)vImageMaxScrCoord.VectorX ;
        for (int jScreen=jStart; jScreen<=jEnd; jScreen++)      {   
            for (int iScreen=iStart; iScreen<=iEnd; iScreen++) {
				this->fImage[jScreen][iScreen] /= this->fMaxImageValue ; 
			}
		}
} 

void Projector::NormalizeRGB(){
        int jStart = (int)vImageMinScrCoord.VectorY ;
        int jEnd   = (int)vImageMaxScrCoord.VectorY ;
        int iStart = (int)vImageMinScrCoord.VectorX ;
        int iEnd   = (int)vImageMaxScrCoord.VectorX ;

		float globalMax= max(max(fMaxRedValue,fMaxGreenValue),fMaxBlueValue) ;
        for (int jScreen=jStart; jScreen<=jEnd; jScreen++)         
			for (int iScreen=iStart; iScreen<=iEnd; iScreen++){
				this->fRGBImage[jScreen][iScreen*4] /= globalMax ;//this->fMaxRedValue ;
				this->fRGBImage[jScreen][iScreen*4+1] /= globalMax ;//this->fMaxGreenValue ;
				this->fRGBImage[jScreen][iScreen*4+2] /= globalMax ; //this->fMaxBlueValue ;
				
			}
			
} 

void Projector::ColorRayCast(Volume& ctvObj)
{
    // set up the FastPass interpolation tables on the first time through
    //if (this->fpiInterpolation == NULL) CalculateFastPassInterpolationTables();

    // set up the pixel step values for fast-pass method vs. normal method
    int jStep                   = 1,    // these are the step values for the normal method.
        iStep                   = 1,    //
        iSimdStep               = 4;    //

    // calculate the two required transformation matrices
    static Transform tObj2Proj,
                     tProj2Obj;
    ctvObj.tCtVolToWorld.LeftMult(this->tWorld2Projector, tObj2Proj);   // make the transform to go from CT space to projector space
    tProj2Obj.Invert(tObj2Proj);                                        // invert it to go from projector space to CT space.
    // map the vertices of the object to projector space and find the limits;
    // these represent the projector frustum coordinates in Projector space
    static Vector vImageMin;
    static Vector vImageMax;
    vImageMin = new Vector( 1e20f, 1e20f, 1e20f);
    vImageMax = new Vector(-1e20f,-1e20f,-1e20f);
    static Vector vObjVertexInCtSpace   (0,0,0),
                  vObjVertexInProjSpace (0,0,0),
                  vObjVertexOnScreen    (0,0,0);
#define FindFrustum(fVertexX,fVertexY,fVertexZ)                                                                      \
    {                                                                                                                \
        vObjVertexInCtSpace.VectorX = fVertexX;                                                                      \
        vObjVertexInCtSpace.VectorY = fVertexY;                                                                      \
        vObjVertexInCtSpace.VectorZ = fVertexZ;                                                                      \
        tObj2Proj.RightMult(vObjVertexInCtSpace,vObjVertexInProjSpace);												 \
		vObjVertexOnScreen.VectorX = vObjVertexInProjSpace.X() * vPrincipalPosition.Z() / vObjVertexInProjSpace.Z(); \
        vObjVertexOnScreen.VectorZ = vPrincipalPosition.Z();														 \
        vObjVertexOnScreen.VectorY = vObjVertexInProjSpace.Y() * vPrincipalPosition.Z() / vObjVertexInProjSpace.Z(); \
        vImageMin.Min(vObjVertexOnScreen);                                                                           \
        vImageMax.Max(vObjVertexOnScreen);                                                                           \
    }
    FindFrustum(ctvObj.vMinNonzeroPlanes.X(), ctvObj.vMinNonzeroPlanes.Y(), ctvObj.vMinNonzeroPlanes.Z());
    FindFrustum(ctvObj.vMaxNonzeroPlanes.X(), ctvObj.vMinNonzeroPlanes.Y(), ctvObj.vMinNonzeroPlanes.Z());
    FindFrustum(ctvObj.vMinNonzeroPlanes.X(), ctvObj.vMaxNonzeroPlanes.Y(), ctvObj.vMinNonzeroPlanes.Z());
    FindFrustum(ctvObj.vMaxNonzeroPlanes.X(), ctvObj.vMaxNonzeroPlanes.Y(), ctvObj.vMinNonzeroPlanes.Z());
    FindFrustum(ctvObj.vMinNonzeroPlanes.X(), ctvObj.vMinNonzeroPlanes.Y(), ctvObj.vMaxNonzeroPlanes.Z());
    FindFrustum(ctvObj.vMaxNonzeroPlanes.X(), ctvObj.vMinNonzeroPlanes.Y(), ctvObj.vMaxNonzeroPlanes.Z());
    FindFrustum(ctvObj.vMinNonzeroPlanes.X(), ctvObj.vMaxNonzeroPlanes.Y(), ctvObj.vMaxNonzeroPlanes.Z());
    FindFrustum(ctvObj.vMaxNonzeroPlanes.X(), ctvObj.vMaxNonzeroPlanes.Y(), ctvObj.vMaxNonzeroPlanes.Z());
#undef FindFrustum
    // all of the above calculations were done in Projector space.  Now we need to calculate the minimum and maximum
    // in screen coordinates and convert to integer values:

    vImageMinScrCoord = new Vector(vImageMin.X() + vCR_Image_Intersection.X(),  // Image X is Projector X, offset by central spot X
                               vImageMin.Y() + vCR_Image_Intersection.Y(),  // Image Y is Projector Z, offset by central spot Y
                               0.0f);                                       // Image Z is irrelevant
    vImageMinScrCoord.Floor();
    vImageMaxScrCoord = new Vector(vImageMax.X() + vCR_Image_Intersection.X(),
                               vImageMax.Y() + vCR_Image_Intersection.Y(),  
                               0.0f);                                       
    vImageMaxScrCoord.Ceil();
    // now we will clip the screen coordinates to the size specified by ipImageSize
    vImageMinScrCoord.VectorX = max(vImageMinScrCoord.VectorX, 0.0f);
    vImageMinScrCoord.VectorY = max(vImageMinScrCoord.VectorY, 0.0f);
    vImageMaxScrCoord.VectorX = min(vImageMaxScrCoord.VectorX, (float)(ipImageSize.x - 1));
    vImageMaxScrCoord.VectorY = min(vImageMaxScrCoord.VectorY, (float)(ipImageSize.y - 1));

    // declare variables to define the frustum limits
    int jStart, jEnd, iStart, iEnd;
    // use one of two raycasting processes, depending on the hardware in THIS machine
   
        // set up accumulators for statistics
        static __m128 m128_Sum, m128_SumSqr, m128_Max;
        m128_Sum        = _mm_setzero_ps();
        m128_SumSqr     = _mm_setzero_ps();
        m128_Max        = _mm_setzero_ps();
		Vector4 v4Max = new Vector4(0.0f, 0.0f, 0.0f) ;			// for maximum calculation

        static int iTupleCnt;
        iTupleCnt       = 0;
        // locate the x-ray source in CT object space
        static Vector  vNull(0,0,0),   // source in projector coordinates; always (0,0,0)
                       vSource;        // source in CT coordinates
        tProj2Obj.RightMult(vNull,vSource);
        // now we can trace each ray through the frustum and add the footprint on the screen
        jStart = (int)vImageMinScrCoord.VectorY,
        jEnd   = (int)vImageMaxScrCoord.VectorY,
        iStart = (int)vImageMinScrCoord.VectorX,
        iEnd   = (int)vImageMaxScrCoord.VectorX;
        for (int jScreen=jStart; jScreen<=jEnd; jScreen+=jStep)         // note that we are stepping wider if bFastPassMethod is 'true'.
            for (int iScreen=iStart; iScreen<=iEnd; iScreen+=iSimdStep) // process four pixels at a whack in the SIMD variables
                                                                        // but skip more X-spaces if bFastPassMethod is 'true'.
        {
            // this is the accumulator that will gather the values of the four pixels
            static __m128 m128_Accum;
			//static __m128 m128_Intensity, m128_Opacity ;
			Vector4 v4Pixel(0.0f, 0.0f, 0.0f) ;

            // first get the locations in Projector space of four X-adjacent pixels
            static Vector4 vPixelInProjectorCoords,
                           vPixelInCtObjCoords;
          
            vPixelInProjectorCoords = new Vector4((0.0f),                                        // this coordinate is filled in below
											  (float)(jScreen - vCR_Image_Intersection.Y()), //Projector Y is image Y  	
											  (float)(vPrincipalPosition.Z()));                // source-to-detector in Projector units
                                             

            for (int i=0; i<4; i++) // load four x-coordinates
            {
                // set the x coordinate of the four pixels
                int ii = iScreen + i * iStep;   // note that we are stepping by two if bFastPassMethod is 'true'
                ((float*)&vPixelInProjectorCoords.VectorX)[i] = (float)(ii - vCR_Image_Intersection.X());// Projector X is minus Image X, 
                                                                                                        // offset by the central spot
				//vPixelInProjectorCoords.VectorX)[i] = (float)(ii - vCR_Image_Intersection.X());// Projector X is minus Image X, 
            }
            // transform the 4 pixels to CT object coordinates
            tProj2Obj.RightMult(vPixelInProjectorCoords,vPixelInCtObjCoords);
            // compute the sampling delta Vector4
            vPixelInCtObjCoords.Subtract(vSource);  // this is the net ray Vector4
            static __m128 m128_RayLength, m128_Normalizer;
            // calculate the magnitudes of the four rays
			// _mm_add_ps  only adds two m128 variables. So, we need to use it twice to add three variables, x, y and z
            m128_RayLength = _mm_sqrt_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(vPixelInCtObjCoords.VectorX,vPixelInCtObjCoords.VectorX),
                                                               _mm_mul_ps(vPixelInCtObjCoords.VectorY,vPixelInCtObjCoords.VectorY)),
                                                               _mm_mul_ps(vPixelInCtObjCoords.VectorZ,vPixelInCtObjCoords.VectorZ) ));
            // calculate a factor to normalize the four rays to unit length (ctvObj.fRayCastStep)
            m128_Normalizer = _mm_div_ps(_mm_set_ps1(ctvObj.fRayCastStep),m128_RayLength);
            // now normalize the four rays to unit length (ctvObj.fRayCastStep)
            vPixelInCtObjCoords.ScalarMult(&m128_Normalizer); // vPixelInCtObjCoords has been converted into the step delta Vector4
#define vRayCastStepDelta   vPixelInCtObjCoords      // we will now use a different name for the same variable
            // Add a 'Kahan zero' to the step delta to avoid a possible divide-by-zero
            vRayCastStepDelta.VectorX = _mm_add_ps(vRayCastStepDelta.VectorX,_mm_set_ps1(1e-20f));
            vRayCastStepDelta.VectorY = _mm_add_ps(vRayCastStepDelta.VectorY,_mm_set_ps1(1e-20f));
            vRayCastStepDelta.VectorZ = _mm_add_ps(vRayCastStepDelta.VectorZ,_mm_set_ps1(1e-20f));
            // find the intersection points of the rays with the bounding planes
            // all of these values are in units of "number of steps"
            static __m128 m128_nx1,m128_nx2,m128_nxL,m128_nxU;
            static __m128 m128_ny1,m128_ny2,m128_nyL,m128_nyU;
            static __m128 m128_nz1,m128_nz2,m128_nzL,m128_nzU;
            {
                m128_nx1 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMinNonzeroPlanes.VectorX),  // x intersection 1
                                                 _mm_set_ps1(vSource.VectorX)                  ),
                                        vRayCastStepDelta.VectorX);
                m128_nx2 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMaxNonzeroPlanes.VectorX),  // x intersection 2
                                                 _mm_set_ps1(vSource.VectorX)                  ),
                                        vRayCastStepDelta.VectorX);
                m128_nxL = _mm_min_ps(m128_nx1,m128_nx2);                                        // first x intersection
                m128_nxU = _mm_max_ps(m128_nx1,m128_nx2);                                        // last x intersection
                m128_ny1 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMinNonzeroPlanes.VectorY),  // y intersection 1
                                                 _mm_set_ps1(vSource.VectorY)                  ),
                                        vRayCastStepDelta.VectorY);
                m128_ny2 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMaxNonzeroPlanes.VectorY),  // y intersection 2
                                                 _mm_set_ps1(vSource.VectorY)                  ),
                                        vRayCastStepDelta.VectorY);
                m128_nyL = _mm_min_ps(m128_ny1,m128_ny2);                                        // first y intersection
                m128_nyU = _mm_max_ps(m128_ny1,m128_ny2);                                        // last y intersection
                m128_nz1 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMinNonzeroPlanes.VectorZ),  // z intersection 1
                                                 _mm_set_ps1(vSource.VectorZ)                  ),
                                        vRayCastStepDelta.VectorZ);
                m128_nz2 = _mm_div_ps(_mm_sub_ps(_mm_set_ps1(ctvObj.vMaxNonzeroPlanes.VectorZ),  // z intersection 2
                                                 _mm_set_ps1(vSource.VectorZ)                  ),
                                        vRayCastStepDelta.VectorZ);
                m128_nzL = _mm_min_ps(m128_nz1,m128_nz2);                                        // first z intersection
                m128_nzU = _mm_max_ps(m128_nz1,m128_nz2);                                        // last z intersection
            }
            // find the start and end points on the four rays, considering the bounding planes
            static __m128 m128_nStart,m128_nEnd;
            static __m128 m128_Count;
            m128_nStart  = _mm_max_ps(m128_nxL,_mm_max_ps(m128_nyL,m128_nzL));    // point of piercing the last of the closest planes (enter)
            m128_nEnd    = _mm_min_ps(m128_nxU,_mm_min_ps(m128_nyU,m128_nzU));    // point of piercing the first of the furthest planes (exit)
            m128_Count   = _mm_sub_ps(m128_nEnd,m128_nStart);                     // number of samples required on each ray
            // if the ray does not pierce the bounded CT object, the count will be <= 0
            int iSampleCount = 0;
            for (int i=0; i<4; i++)
            {
                // calculate the maximum number of samples we must take in order to get the longest ray
                iSampleCount = max(iSampleCount, (int)((float*)&m128_Count)[i]);  // drop fractional part
            }
            // calculate the first sample point for each ray
            static Vector4 vSamplePoint;
            vSamplePoint = vRayCastStepDelta;
            if (iSampleCount <= 0) goto NoPiercing; // we can skip these four rays completely if none pierce the CT object!
            vSamplePoint.ScalarMult(&m128_nStart);
            vSamplePoint.Add(vSource);
            // now we can run the rays
            for (int iSample=0; iSample<=iSampleCount; iSample++)
            {
                // round the indices to nearest neighbor
                static __m128i m128i_i,m128i_j,m128i_k;
                m128i_i = _mm_cvtps_epi32(vSamplePoint.VectorX);
                m128i_j = _mm_cvtps_epi32(vSamplePoint.VectorY);
                m128i_k = _mm_cvtps_epi32(vSamplePoint.VectorZ);
                // accumulate an __m128 of voxels by picking the four values out of the CT object
                for (int i=0; i<4; i++)
                {
                    if (((float*)&m128_Count)[i] >= iSample)
                    {
                        assert(((int*)&m128i_k)[i] >= ctvObj.vMinNonzeroPlanes.VectorZ && 
                               ((int*)&m128i_k)[i] <= ctvObj.vMaxNonzeroPlanes.VectorZ);
                        assert(((int*)&m128i_j)[i] >= ctvObj.vMinNonzeroPlanes.VectorY && 
                               ((int*)&m128i_j)[i] <= ctvObj.vMaxNonzeroPlanes.VectorY);
                        assert(((int*)&m128i_i)[i] >= ctvObj.vMinNonzeroPlanes.VectorX && 
                               ((int*)&m128i_i)[i] <= ctvObj.vMaxNonzeroPlanes.VectorX);

						//m128_Intensity)[i] =  (float)ctvObj.usCtVol[m128i_k)[i]]
                        //                                                [m128i_j)[i]]
                        //                                                [m128i_i)[i]] ;

						float intensity = (float)ctvObj.usCtVol[((int*)&m128i_k)[i]]
                                                                        [((int*)&m128i_j)[i]]
                                                                        [((int*)&m128i_i)[i]] ;

																		
						float opacity = this->opacityFunction.getValue( intensity );

						//printf("intensity: %f opacity: %f\n", intensity, opacity);

						__m128 pixelRGB = this->colorTranFun.getValue(intensity, opacity) ;  
						((float*)&v4Pixel.VectorX)[i] += ((float*)&pixelRGB)[0] ;		// RED
						((float*)&v4Pixel.VectorY)[i] += ((float*)&pixelRGB)[1] ;		// GREEN
						((float*)&v4Pixel.VectorZ)[i] += ((float*)&pixelRGB)[2] ;		// BLUE

                    }					
                }

                // index to the next four sample positions
                vSamplePoint.VectorX = _mm_add_ps(vSamplePoint.VectorX,vRayCastStepDelta.VectorX);
                vSamplePoint.VectorY = _mm_add_ps(vSamplePoint.VectorY,vRayCastStepDelta.VectorY);
                vSamplePoint.VectorZ = _mm_add_ps(vSamplePoint.VectorZ,vRayCastStepDelta.VectorZ);
            }

            // put away the rays differently, depending on bFastPassMethod switch
            /*if (bFastPassMethod)
            {
                fImage[jScreen][iScreen                                ] = m128_Accum)[0];
                fImage[jScreen][iScreen +     iFastPassSamplingInterval] = m128_Accum)[1];
                fImage[jScreen][iScreen + 2 * iFastPassSamplingInterval] = m128_Accum)[2];
                fImage[jScreen][iScreen + 3 * iFastPassSamplingInterval] = m128_Accum)[3];
            } */  // if (this->bFastPassMethod)
            //else    // if ( ! this->bFastPassMethod) //i.e., normal method
            {
				int iEffScreen = iScreen*4 ;
                fRGBImage[jScreen][iEffScreen + 0] = ((float*)&v4Pixel.VectorX)[0];		//RED 
				fRGBImage[jScreen][iEffScreen + 1] = ((float*)&v4Pixel.VectorY)[0];		//GREEN 
				fRGBImage[jScreen][iEffScreen + 2] = ((float*)&v4Pixel.VectorZ)[0];		//BLUE 
				fRGBImage[jScreen][iEffScreen + 3] = 0.5;								//ALPHA
				//second pixel
                fRGBImage[jScreen][iEffScreen +  4] = ((float*)&v4Pixel.VectorX)[1];		//RED 
				fRGBImage[jScreen][iEffScreen +  5] = ((float*)&v4Pixel.VectorY)[1];		//GREEN 
				fRGBImage[jScreen][iEffScreen +  6] = ((float*)&v4Pixel.VectorZ)[1];		//BLUE
				fRGBImage[jScreen][iEffScreen +  7] = 0.5;								//ALPHA
				// third pixel
                //fRGBImage[jScreen][iScreen+2    ][0] = v4Pixel.VectorX)[2];		//RED 
				//fRGBImage[jScreen][iScreen+2    ][1] = v4Pixel.VectorY)[2];		//GREEN 
				//fRGBImage[jScreen][iScreen+2    ][2] = v4Pixel.VectorZ)[2];		//BLUE
				//fRGBImage[jScreen][iScreen+2    ][3] = 0.5;								//ALPHA

				// third pixel
                fRGBImage[jScreen][iEffScreen +  8] = ((float*)&v4Pixel.VectorX)[2];		//RED 
				fRGBImage[jScreen][iEffScreen +  9] = ((float*)&v4Pixel.VectorY)[2];		//GREEN 
				fRGBImage[jScreen][iEffScreen + 10] = ((float*)&v4Pixel.VectorZ)[2];		//BLUE
				fRGBImage[jScreen][iEffScreen + 11] = 0.5;								//ALPHA

				// fourth pixel
                fRGBImage[jScreen][iEffScreen + 12] = ((float*)&v4Pixel.VectorX)[3];		//RED 
				fRGBImage[jScreen][iEffScreen + 13] = ((float*)&v4Pixel.VectorY)[3];		//GREEN 
				fRGBImage[jScreen][iEffScreen + 14] = ((float*)&v4Pixel.VectorZ)[3];		//BLUE
				fRGBImage[jScreen][iEffScreen + 15] = 0.5;								//ALPHA

            }   // else    // if ( ! this->bFastPassMethod) //i.e., normal method

            // accumulate the stats
            //m128_Sum    = _mm_add_ps(m128_Sum,m128_Accum);
            //m128_SumSqr = _mm_add_ps(m128_SumSqr,_mm_mul_ps(m128_Accum,m128_Accum));

			v4Max.VectorX = _mm_max_ps(v4Max.VectorX , v4Pixel.VectorX);		// red
			v4Max.VectorY = _mm_max_ps(v4Max.VectorY , v4Pixel.VectorY);		// green
			v4Max.VectorZ = _mm_max_ps(v4Max.VectorZ , v4Pixel.VectorZ);		// blue
            iTupleCnt++;
NoPiercing:
            continue;
        }   // for (int jScreen=jStart; jScreen<=jEnd; jScreen+=jStep) for (int iScreen=iStart; iScreen<=iEnd; iScreen+=iSimdStep)
        // finish the stats
        //fMeanImageValue          = (float)(AccumulateFour(m128_Sum)    / (4.0f * iTupleCnt));
        //float fMeanSqrImageValue = (float)(AccumulateFour(m128_SumSqr) / (4.0f * iTupleCnt));
        //fStdDevImageValue        = sqrtf(fMeanSqrImageValue - sqr(fMeanImageValue));
		this->fMaxRedValue  = MaxFour(v4Max.VectorX) ;
		this->fMaxGreenValue = MaxFour(v4Max.VectorY) ;
		this->fMaxBlueValue = MaxFour(v4Max.VectorZ) ;

#ifdef WriteProjectionFile
    {
        char szFileName[MAX_PATH];
        sprintf(szFileName,
                "%sDRRbeforeInterpolation.raw",
                this->szProjectorName);
        FILE* fpImageFile = fopen(szFileName, "wb");
        if (fpImageFile)
        {
            fwrite( &fImage[0][0], 
                    sizeof(float), 
                    ipImageSize.x * ipImageSize.y, 
                    fpImageFile);
            fclose( fpImageFile);
        }
        else
        {
            FATAL_ERROR(szFileName);
        }
    }   // #ifdef WriteProjectionFile
#endif 

#ifdef WriteProjectionFile
    {
        char szFileName[MAX_PATH];
        sprintf(szFileName,
                "%sDRRafterInterpolation.raw",
                this->szProjectorName);
        FILE* fpImageFile = fopen(szFileName, "wb");
        if (fpImageFile)
        {
            fwrite( &fImage[0][0], 
                    sizeof(float), 
                    ipImageSize.x * ipImageSize.y, 
                    fpImageFile);
            fclose( fpImageFile);
        }
        else
        {
            FATAL_ERROR(szFileName);
        }
    }   // #ifdef WriteProjectionFile
#endif 

    // mark the end time and accumulate it

	// now we will extend the boundary if there was any more projection already in the frameBuffer
//    DWORD dwEndTime = GetTickCount();
//  this->iLastProjectionTime   = (dwEndTime - dwStartTime);
//    this->iTotalProjectionTime += this->iLastProjectionTime;
    this->iProjectionCount++;
}


// Return the average number of milliseconds used per rendering operation
int Projector::AverageRenderTime(void)
{
    if (this->iProjectionCount)
    {
        return this->iTotalProjectionTime / this->iProjectionCount;
    }
    else
    {
        return 0;
    }
}

// Calculates the world-to-projector and projector-to-world transforms from the position and orientation members
void Projector::CalculateTransformations(void)
{
    this->vPrincipalPosition        = new Vector(0.0f,
                                             0.0f,
                                             -(this->fSourceToDetector / this->fDetectorPixelSize) ); // since projector projects towards negative z axis

    Transform tScale,tMove,tXrot,tYrot,tZrot;                   // set up the individual transformations
    tScale.Scale(vProjectorScale);
    tZrot.RotateZ(vCentralRayOrientation.Z() * M_PI / 180.0);
    tYrot.RotateY(vCentralRayOrientation.Y() * M_PI / 180.0);
    tXrot.RotateX(vCentralRayOrientation.X() * M_PI / 180.0);
    tMove.Translate(this->vSourceLocation );
    tProjector2World = tScale;                                  // apply the transforms in sequence
    tProjector2World.LeftMult(tZrot);
    tProjector2World.LeftMult(tYrot);
    tProjector2World.LeftMult(tXrot);
    tProjector2World.LeftMult(tMove);
    tWorld2Projector.Invert(tProjector2World);                  // invert the transform for the reverse map
}

void Projector::CalculateFastPassInterpolationTables(void)
{
    // allocate the 3D array for fast-pass interpolation
    this->fpiInterpolation = (FPinterp***)AllocArray3D( 5,                                      // interpolation elements present
                                                        this->iFastPassSamplingInterval + 1,    // j index being interpolated
                                                        this->iFastPassSamplingInterval + 1,    // i index being interpolated
                                                        sizeof(FPinterp));                      // fpWeight[jj][ii]

    // fpiInterpolation [k][j][i]:
    // k=0 ==> Sample[0][0] is not valid; i.e. RayCast() did not project a value through the CT.
    // k=1 ==> Sample[0][1] is not valid; i.e. RayCast() did not project a value through the CT.
    // k=2 ==> Sample[1][0] is not valid; i.e. RayCast() did not project a value through the CT.
    // k=3 ==> Sample[1][1] is not valid; i.e. RayCast() did not project a value through the CT.
    // k=4 ==> All four corners are valid samples; this is standard bi-linear interpolation.

    // first calculate the standard bi-linear weights
    for (int j=0; j<=iFastPassSamplingInterval; j++) 
        for (int i=0; i<=iFastPassSamplingInterval; i++) 
    {
        double fJweight = (double)j / (double)iFastPassSamplingInterval;
        double fIweight = (double)i / (double)iFastPassSamplingInterval;
        double dSum     = 0;
        dSum += fpiInterpolation[4][j][i].fpWeight[0][0] =  (float)((1.0 - fJweight) * (1.0 - fIweight));
        dSum += fpiInterpolation[4][j][i].fpWeight[0][1] =  (float)((1.0 - fJweight) * (      fIweight));
        dSum += fpiInterpolation[4][j][i].fpWeight[1][0] =  (float)((      fJweight) * (1.0 - fIweight));
        dSum += fpiInterpolation[4][j][i].fpWeight[1][1] =  (float)((      fJweight) * (      fIweight));

        // the four weights should add to 1 within single-precision limits.
        if (fabs(dSum - 1.0) > 0.000001) FATAL_ERROR("Unable to build the FastPass bi-linear table");

    }   // for (int j=0; j<=iFastPassSamplingInterval; j++) for (int i=0; i<=iFastPassSamplingInterval; i++) 

    // calculate the partial cases:
    for (int j=0; j<=iFastPassSamplingInterval; j++) 
        for (int i=0; i<=iFastPassSamplingInterval; i++) 
    {
        // compute k==0 (Sample[0][0] missing)
        if ( (i + j) >= iFastPassSamplingInterval)
        {
            // we need to adjust out the contribution of Sample[0][0]
            double dSm4 =   // fpiInterpolation[4][j][i].fpWeight[0][0] +
                            fpiInterpolation[4][j][i].fpWeight[0][1] +
                            fpiInterpolation[4][j][i].fpWeight[1][0] +
                            fpiInterpolation[4][j][i].fpWeight[1][1] ;
            double dSum = 0;
            dSum += fpiInterpolation[0][j][i].fpWeight[0][0] = 0;   // (float)(fpiInterpolation[4][j][i].fpWeight[0][0] / dSm4);
            dSum += fpiInterpolation[0][j][i].fpWeight[0][1] = (float)(fpiInterpolation[4][j][i].fpWeight[0][1] / dSm4);
            dSum += fpiInterpolation[0][j][i].fpWeight[1][0] = (float)(fpiInterpolation[4][j][i].fpWeight[1][0] / dSm4);
            dSum += fpiInterpolation[0][j][i].fpWeight[1][1] = (float)(fpiInterpolation[4][j][i].fpWeight[1][1] / dSm4);

            // the four weights should add to 1 within single-precision limits.
            if (fabs(dSum - 1.0) > 0.000001) FATAL_ERROR("Unable to build the FastPass k==0 table");

        }   // if ( (i + j) >= iFastPassSamplingInterval)

        // compute k==1 (Sample[0][1] missing)
        if ( (iFastPassSamplingInterval - i + j) >= iFastPassSamplingInterval)
        {
            // we need to adjust out the contribution of Sample[0][1]
            double dSm4 =   fpiInterpolation[4][j][i].fpWeight[0][0] +
                            // fpiInterpolation[4][j][i].fpWeight[0][1] +
                            fpiInterpolation[4][j][i].fpWeight[1][0] +
                            fpiInterpolation[4][j][i].fpWeight[1][1] ;
            double dSum = 0;
            dSum += fpiInterpolation[1][j][i].fpWeight[0][0] = (float)(fpiInterpolation[4][j][i].fpWeight[0][0] / dSm4);
            dSum += fpiInterpolation[1][j][i].fpWeight[0][1] = 0; // (float)(fpiInterpolation[4][j][i].fpWeight[0][1] / dSm4);
            dSum += fpiInterpolation[1][j][i].fpWeight[1][0] = (float)(fpiInterpolation[4][j][i].fpWeight[1][0] / dSm4);
            dSum += fpiInterpolation[1][j][i].fpWeight[1][1] = (float)(fpiInterpolation[4][j][i].fpWeight[1][1] / dSm4);

            // the four weights should add to 1 within single-precision limits.
            if (fabs(dSum - 1.0) > 0.000001) FATAL_ERROR("Unable to build the FastPass k==1 table");

        }   // if ( (iFastPassSamplingInterval - i + j) >= iFastPassSamplingInterval)

        // compute k==2 (Sample[1][0] missing)
        if ( (i + iFastPassSamplingInterval - j) >= iFastPassSamplingInterval)
        {
            // we need to adjust out the contribution of Sample[1][0]
            double dSm4 =   fpiInterpolation[4][j][i].fpWeight[0][0] +
                            fpiInterpolation[4][j][i].fpWeight[0][1] +
                            // fpiInterpolation[4][j][i].fpWeight[1][0] +
                            fpiInterpolation[4][j][i].fpWeight[1][1] ;
            double dSum = 0;
            dSum += fpiInterpolation[2][j][i].fpWeight[0][0] = (float)(fpiInterpolation[4][j][i].fpWeight[0][0] / dSm4);
            dSum += fpiInterpolation[2][j][i].fpWeight[0][1] = (float)(fpiInterpolation[4][j][i].fpWeight[0][1] / dSm4);
            dSum += fpiInterpolation[2][j][i].fpWeight[1][0] = 0; // (float)(fpiInterpolation[4][j][i].fpWeight[1][0] / dSm4);
            dSum += fpiInterpolation[2][j][i].fpWeight[1][1] = (float)(fpiInterpolation[4][j][i].fpWeight[1][1] / dSm4);

            // the four weights should add to 1 within single-precision limits.
            if (fabs(dSum - 1.0) > 0.000001) FATAL_ERROR("Unable to build the FastPass k==2 table");

        }   // if ( (i + iFastPassSamplingInterval - j) >= iFastPassSamplingInterval)

        // compute k==3 (Sample[1][1] missing)
        if ( (iFastPassSamplingInterval * 2 - i - j) >= iFastPassSamplingInterval)
        {
            // we need to adjust out the contribution of Sample[1][1]
            double dSm4 =   fpiInterpolation[4][j][i].fpWeight[0][0] +
                            fpiInterpolation[4][j][i].fpWeight[0][1] +
                            fpiInterpolation[4][j][i].fpWeight[1][0] ;
                            // fpiInterpolation[4][j][i].fpWeight[1][1] ;
            double dSum = 0;
            dSum += fpiInterpolation[3][j][i].fpWeight[0][0] = (float)(fpiInterpolation[4][j][i].fpWeight[0][0] / dSm4);
            dSum += fpiInterpolation[3][j][i].fpWeight[0][1] = (float)(fpiInterpolation[4][j][i].fpWeight[0][1] / dSm4);
            dSum += fpiInterpolation[3][j][i].fpWeight[1][0] = (float)(fpiInterpolation[4][j][i].fpWeight[1][0] / dSm4);
            dSum += fpiInterpolation[3][j][i].fpWeight[1][1] = 0; // (float)(fpiInterpolation[4][j][i].fpWeight[1][1] / dSm4);

            // the four weights should add to 1 within single-precision limits.
            if (fabs(dSum - 1.0) > 0.000001) FATAL_ERROR("Unable to build the FastPass k==3 table");

        }   // if ( (iFastPassSamplingInterval * 2 - i - j) >= iFastPassSamplingInterval)

    }   // for (int j=0; j<=iFastPassSamplingInterval; j++) for (int i=0; i<=iFastPassSamplingInterval; i++) 

#ifdef DumpFastPassInterpolationFactors
    FILE* fpDump = fopen("FastPassInterpolationFactors.csv", "wt");
    if (fpDump == NULL) FATAL_ERROR("FastPassInterpolationFactors.csv");

    for (int k=0; k<5; k++) 
        for (int jj=0; jj<2; jj++) 
            for (int ii=0; ii<2; ii++) 
    {
        fprintf(fpDump, "\nk=%d,jj=%d,ii=%d\n", k, jj, ii);
        for (int j=0; j<=iFastPassSamplingInterval; j++)
        {
            for (int i=0; i<=iFastPassSamplingInterval; i++)
            {
                fprintf(fpDump, "%f,", fpiInterpolation[k][j][i].fpWeight[jj][ii]);
            }   // for (int i=0; i<=iFastPassSamplingInterval; i++)
            fprintf(fpDump,"\n");
        }   // for (int j=0; j<=iFastPassSamplingInterval; j++)
    }   // for (int k=0; k<5; k++) for (int jj=0; jj<2; jj++) for (int ii=0; ii<2; ii++) 

    fclose(fpDump);

#endif  // #ifdef DumpFastPassInterpolationFactors
}

void Projector::printProjector(void) {
	printf("Projector Name: %s\n",this->szProjectorName) ;
	printf("Projector location %s\n",this->vSourceLocation.Format()) ;
	printf("Projector central ray orientation %s\n",this->vCentralRayOrientation.Format() ) ;
	printf("Projector VCR_Image_Intersection %s\n",this->vCR_Image_Intersection.Format()  ) ;
	printf("Projector scale %s\n",this->vProjectorScale.Format()  ) ;
	printf("Source to Detector %f\n",this->fSourceToDetector) ;
	printf("Detector pixel size %f\n", this->fDetectorPixelSize) ;
	printf("Projector principal position %s\n",this->vPrincipalPosition.Format()  ) ;
	printf("Projector min screen coordinate %s\n",this->vImageMinScrCoord.Format()  ) ;
	printf("Projector max screen coordinate %s\n",this->vImageMaxScrCoord.Format()  ) ;
	printf("Transform matrix for world to projector\n") ;
	this->tWorld2Projector.Print() ;
	printf("Bead count %d\n",this->iBeadCount) ;
	printf("Calibaration error %f\n",this->dCalibrationError ) ;	
}
