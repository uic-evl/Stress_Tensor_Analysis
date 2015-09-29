/* Taken from MTW project by Abedul Haque, December 2009		*/
/* This file has been modified almost 90%						*/
/* The original file was much more bigger and complicated		*/
/* I Kept what was necessary for volume rendering settings.		*/
/* Comments may not seem appropriate since I didn't change 'em	*/
/* used for CS3610 class project								*/

#pragma once

#include "../vector/Point.hh"
#include "../vector/Vector.hh"
#include "../transform/TransferFunction.hh"
#include "../volume/ColorTransferFunction.hh"
#include <iostream>

#define ParameterSeparators                 " ,\t/\n;"
class Parameters
{
// member variables
public:
// Volume file parameters
    char    CtFileName[MAX_PATH];               			// Volume file to be matched
	char	szFilePath[MAX_PATH];
    iPoint  CtFileSize;                         // Volume file dimensions
    Vector  CtFileVoxelSpacing;                 // voxel spacing of file
    Vector  CtStartingPosition;                 // allows the operator to specify a starting position for the CT object at load time
    Vector  CtStartingOrientation;              // allows the operator to specify a starting orientation for the CT object at load time
    float   CtSamplingInterval;                 // distance between ray-cast samples in millimeters

// red fluoroscope parameters
    Vector  RedSourceLoc;                       // location of x-ray source in world coordinates (mm)
    Vector  RedOrientation;                     // central ray orientation vs. world axes (just a made-up configuration)
    float   RedSourceToDetector;                // source-to-detector distance (in mm)
    float   RedPixelSize;                       // detector pixel size (in mm)
    Vector  RedCentralSpotOnScreen;             // location of central spot in screen coordinates
    iPoint  RedMovieFileSize;                   // dimensions of the fluoroscope file
                                                // NOTE: additional codes have been created (05-24-2007)

	TransferFunction opTranFunc;				// Opacity transfer function of the projector
	ColorTransferFunction colTranFunc ;			// color transfer function for the projector

	fPoint  SeedPointLoc;

#define MTW_Default_FastPassSamplingInterval    8
    int     FastPassSamplingInterval;           // FastPass mode samples the CT DRR at this interval in both X and Y

// member functions
public:
    Parameters(void);
    ~Parameters(void);
    // Read in all operational parameters from the Parameters.csv file
    void LoadParameters(std::string szpFile);
    // process BOTH of the column-order variables
    void ProcessColumnOrderCommand(void);
    // Write out all operational parameters in a new Parameters*.csv file
    void SaveParameters(void);
};
