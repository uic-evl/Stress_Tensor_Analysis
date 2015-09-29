/* this file has been taken from MTW project		*/
/* I modified the file as needed. There have been	*/
/* a lot of changes. I didn't change the comments.	*/
/* So, some comments may be misleading.				*/		


// This class supports the loading and handling of a volume of CT-scan data.
// The data are assumed to be 16-bit, unsigned integers.  The volumes are
// of various dimensions, and are presented without header information (raw binary).
// The voxel spacing is assumed to be the same in the X and Y directions, but
// usually different in the Z direction (slice spacing).

#define _USE_MATH_DEFINES

#include <iostream>
#include "math.h"
#include "Volume.hh"
#include "../std/AllocArray.hh"
#include "../std/FatalError.hh"

Volume::Volume(void)
{
    memset(this,0,sizeof(Volume));
}

Volume::~Volume(void)
{}


void Volume::ComputeTransforms(void)
{
    // calculate the individual moves and rotates
    Transform tObj2CtVol,tScale,tXrot,tYrot,tZrot,tPosition;
    Vector vCtVol000 = vOrigin;
    vCtVol000.Negate();
    tObj2CtVol.Translate(vCtVol000);

    tScale.Scale(vVoxelSpacing);
    tXrot.RotateX(vOrientation.X() * (float)(M_PI / 180.0));
    tYrot.RotateY(vOrientation.Y() * (float)(M_PI / 180.0));
    tZrot.RotateZ(vOrientation.Z() * (float)(M_PI / 180.0));
    tPosition.Translate(vPosition);
    // now stack them up to build one transform matrix
    tCtVolToWorld = tObj2CtVol;
    tCtVolToWorld.LeftMult(tScale);
    tCtVolToWorld.LeftMult(tZrot);
    tCtVolToWorld.LeftMult(tYrot);
    tCtVolToWorld.LeftMult(tXrot);
    tCtVolToWorld.LeftMult(tPosition);
    // now invert the matrix to get the world-to-CT transform
    tWorldToCtVol.Invert(tCtVolToWorld);
#ifdef _DEBUG
    printf("\nTest transformations from CT volume to world coordinates:\n");
    Vector vCT(vOrigin),vWorld;
    tCtVolToWorld.RightMult(vCT,vWorld);
    vCT.Print();
    vWorld.Print();
    vCT.VectorX += 1.0;
    tCtVolToWorld.RightMult(vCT,vWorld);
    vCT.Print();
    vWorld.Print();
    vCT.VectorY += 1.0;
    tCtVolToWorld.RightMult(vCT,vWorld);
    vCT.Print();
    vWorld.Print();
    vCT.VectorZ += 1.0;
    tCtVolToWorld.RightMult(vCT,vWorld);
    vCT.Print();
    vWorld.Print();
#endif // _DEBUG
}

// Sets the rotational origin of the volume to vNewOrigin (in file coordinates)
void Volume::SetVolumeOrigin(Vector& vNewOrigin)
{
    // load the new file origin
    vFileOrigin = vNewOrigin;
    // calculate the interpolated volume origin
    vOrigin = new Vector((float)(vFileOrigin.X() * ipStretch.x),
                     (float)(vFileOrigin.Y() * ipStretch.y),
                     (float)(vFileOrigin.Z() * ipStretch.z));
    // report the chante
    printf("Rotational origin of CT file changed to %s file voxels.\n",vFileOrigin.Format());
}

// Load the CT volume from a source file without any stretching.
void Volume::LoadCtSansInterpolation(char*   szInFileName, 
                                       iPoint& ipDimensions, 
                                       Vector& vVoxelSpacingOnFile, 
                                       float   fCtSamplingInterval)
{
    int i,j,k;     // utility indices
    // copy the parameters into this object
    strcpy(this->szFileName,szInFileName);
    this->ipFileSize        = ipDimensions;
    this->vFileVoxelSpacing = vVoxelSpacingOnFile;
    this->ipStretch         = iPoint(1,1,1);            // these data will be ignored anyway.
    // allocate the main volume image array
    this->ipVolSize = ipDimensions; // the volume dimensions are the same as the input CT file
    this->usCtVol   = /*new u_short**[ipVolSize.x];*/(u_short***)AllocArray3D(ipVolSize.z,
                                               ipVolSize.y,
                                               ipVolSize.x,
                                               sizeof(u_short)); 
/*
	for(int ii = 0; ii <= ipVolSize.x; ii++){
		this->usCtVol[ii] = new u_short*[ipVolSize.y];
		for(int jj = 0; jj <= ipVolSize.y; jj++)
			this->usCtVol[ii][jj] = new u_short[ipVolSize.z];
	}
*/

    {
        // open the input file
        FILE* fpIn = fopen(szFileName,"rb");
        if (fpIn)
        {
            printf("'%s' opened for input\n",szFileName);
        }
        else
        {
            FATAL_ERROR("std error\n");
        }

        // for(i = 0; i < ipVolSize.x; i++)
        // {
        //     for (j = 0; j < ipVolSize.y; j++)
        //     {
        //         for (k = 0; k< ipVolSize.z; k++)
        //         {      
        //             u_short temp;
        //                         size_t stRead = fread(  &temp,
        //                             sizeof(u_short),
        //                             1,
        //                             fpIn);
        //             std::cout << temp << std::endl;
        //         }
        //     }
        // }

        // read the file one plane at a time and load them in reverse order so as to match the projector geometry
        for (k=0; k<=(ipVolSize.z-1); k++)
        {
            size_t stRead = fread(  &usCtVol[k][0][0],
                                    sizeof(u_short),
                                    ipVolSize.x * ipVolSize.y,
                                    fpIn);
            std::cout << stRead << std::endl;
            
            if ((int)stRead != (ipVolSize.x * ipVolSize.y))
            {
                printf("%d\n", (int)ipVolSize.x * ipVolSize.y );
                // short read -- this is a problem
                FATAL_ERROR("Premature end-of-file on the CT volume");
            }
        }

        fclose(fpIn);
    }   // else// if (! bFileIsTiff) 

	vOrigin = new Vector( (float)(ipVolSize.x)/2.0, 
						(float)(ipVolSize.y)/2.0,
						(float)(ipVolSize.z)/2.0) ;
    vMinNonzeroPlanes = new Vector(0.0,0.0,0.0);
    vMaxNonzeroPlanes = new Vector(ipVolSize.x-1 ,ipVolSize.y-1 ,ipVolSize.z-1);

    // report the origin in terms of the original file voxels
    vFileOrigin = new Vector(vOrigin.X(),
                         vOrigin.Y(),
                         ipFileSize.z - vOrigin.Z());
    //printf("Rotational origin of CT file is @ %s in file voxels.\n",vFileOrigin.Format());
    // calculate the voxel spacing in the stored object
    vVoxelSpacing = new Vector(vFileVoxelSpacing.X() / ipStretch.x,
                           vFileVoxelSpacing.Y() / ipStretch.y,
                           vFileVoxelSpacing.Z() / ipStretch.z);
    // calculate the default ray-casting step interval in units of the usCtVol[][][] array
    if (fCtSamplingInterval)
    {
        this->fRayCastStep = fCtSamplingInterval / max(max(vVoxelSpacing.X(),
                                                           vVoxelSpacing.Y()),
                                                           vVoxelSpacing.Z());
    }
    else
    {
        this->fRayCastStep = (float)min(ipStretch.x, min(ipStretch.y, ipStretch.x));
    }
    // clear the position and orientation to zero
    vPosition = vOrientation = new Vector(0,0,0);
 
    // mark the volume as loaded
    bCtVolumeIsLoaded = true;
}
