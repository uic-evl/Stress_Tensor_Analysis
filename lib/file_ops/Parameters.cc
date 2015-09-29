/* Taken from MTW project by Abedul Haque, December 2009		*/
/* This file has been modified almost 90%						*/
/* The original file was much more bigger and complicated		*/
/* I Kept what was necessary for volume rendering settings.		*/
/* Comments may not seem appropriate since I didn't change 'em	*/
/* used for CS3610 class project								*/

#include "Parameters.hh"
#include "../std/FatalError.hh"
#include <cstdio>
#include <cstdlib>
#include <string.h>

extern int slices;
extern int rows;
extern int cols;

Parameters::Parameters(void)
{
    // wipe the control block to all zeros
    memset(this,0,sizeof(Parameters));
    // initialize default values
    int i,j;
    CtFileVoxelSpacing      = new Vector(0,0,0);    // initialize the standard Vector objects to (0,0,0)
    CtStartingPosition      = new Vector(0,0,0);
    CtStartingOrientation   = new Vector(0,0,0);

    RedSourceLoc            = new Vector(0,0,0);
    RedOrientation          = new Vector(0,0,0);
    RedCentralSpotOnScreen  = new Vector(0,0,0);

    FastPassSamplingInterval= MTW_Default_FastPassSamplingInterval;
}

Parameters::~Parameters(void)
{
}

// Read in all operational parameters from the Parameters.csv file
void Parameters::LoadParameters(std::string szpFile)
{
    char szWork[MAX_PATH * 10]; // utility string buffer

    // open the parameter file
//    strcpy(this->szFilePath, szpFile);   // file name from the command line
    FILE* fpParm = fopen(szpFile.c_str(),"rt");
//    if (strlen(szFilePath) > 0) fpParm = fopen(szFilePath,"rt");
    char szLine[MAX_PATH];
    if (fpParm == NULL)
    {
        // we couldn't find the file -- ask the operator for a file path using the common controls file dialog
		FATAL_ERROR("Unable to continue without a parameter file.");
    }
    // read and process each line
#define IsKeyWord(x) (strcmp(szpKeyWord,x) == 0)
    for (int iLine=2; iLine<1000; iLine++)
    {
        if (fgets(szLine,MAX_PATH,fpParm) == NULL )  break;  // end of file
        char* szpKeyWord = strtok(szLine,ParameterSeparators);
        if (szpKeyWord == NULL) continue;   // skip blank lines
//        _strupr(szpKeyWord);
        // isolate the parameters on the line
        char* szpParm1 = strtok(NULL,ParameterSeparators);
        char* szpParm2 = strtok(NULL,ParameterSeparators);
        char* szpParm3 = strtok(NULL,ParameterSeparators);
		char* szpParm4 = strtok(NULL,ParameterSeparators);

        if (IsKeyWord("KeyWord") || IsKeyWord("Header"))
        {   
            continue;   // this line is just a comment
        }
        if (IsKeyWord("CtFileName"))
        {   	
			std:: cout << szpFile << std::endl;
			int pos = szpFile.find("parameters.csv");
			std::string path = szpFile.substr(0, pos);
			
			path.append("divergence.raw");
			
			strncpy(CtFileName, path.c_str(), path.size());
						
        }
        else if (IsKeyWord("CtFileSize"))
        {   
                CtFileSize = iPoint(atoi(szpParm1),atoi(szpParm2),atoi(szpParm3));
				cols = atoi(szpParm1);
				rows = atoi(szpParm2);
				slices = atoi(szpParm3);
        }
        else if (IsKeyWord("CtFileVoxelSpacing"))
        {   
            CtFileVoxelSpacing = new Vector(atof(szpParm1),atof(szpParm2),atof(szpParm3));
        }
        else if (IsKeyWord("CtSamplingInterval"))
        {   
            CtSamplingInterval = (float)atof(szpParm1);
        }
        else if (IsKeyWord("RedSourceLoc"))
        {   
            RedSourceLoc = new Vector(atof(szpParm1),atof(szpParm2),atof(szpParm3));
        }
        else if (IsKeyWord("RedOrientation"))
        {   
            RedOrientation = new Vector(atof(szpParm1),atof(szpParm2),atof(szpParm3));
        }
        else if (IsKeyWord("RedSourceToDetector"))
        {   
            RedSourceToDetector = (float)atof(szpParm1);
        }
        else if (IsKeyWord("RedPixelSize"))
        {   
            // we only accept the pixel size parameter if it was not encoded in the movie file
            if (RedPixelSize == 0) RedPixelSize = (float)atof(szpParm1);
        }
        else if (IsKeyWord("RedCentralSpotOnScreen"))
        {   
            RedCentralSpotOnScreen = new Vector(atof(szpParm1),atof(szpParm2),0.0);
        }
        else if (IsKeyWord("RedMovieFileSize"))
        {   
                RedMovieFileSize = iPoint(atoi(szpParm1),atoi(szpParm2),atoi(szpParm3));
        }
        else if (IsKeyWord("FastPassSamplingInterval"))
        {   
            FastPassSamplingInterval = atoi(szpParm1);
        }
		else if(IsKeyWord("IntensityTranFunc")){
			this->opTranFunc.addPoint(atof(szpParm1), atof(szpParm2)) ;
		}
		else if(IsKeyWord("ColorTranFunc")){
			this->colTranFunc.addPoint(atof(szpParm1), atof(szpParm2), atof(szpParm3), atof(szpParm4)) ;
		}
        else if (IsKeyWord("SeedPoint"))
        {   
			this->SeedPointLoc = fPoint(atof(szpParm1),atof(szpParm2),atof(szpParm3));
			//printf("%f %f %f\n", SeedPointLoc.x,SeedPointLoc.y,SeedPointLoc.z );
        }
        else
        {
            sprintf(szWork,
                    "Unable to interpret keyword '%s' in line %d of the parameter file\n%s",
                    szpKeyWord,
                    iLine,
                   szpFile.c_str());
            FATAL_ERROR(szWork); 
        }
    }
    // close the file and quit
    fclose(fpParm);

    // check the pixel size
//    if (this->RedPixelSize   == 0)  FATAL_ERROR("You must enter the pixel size for the incoming fluoroscope files.");
 }

