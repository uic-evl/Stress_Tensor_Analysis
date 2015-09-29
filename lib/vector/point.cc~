// point.cpp

// routines to support iPoint, fPoint and dPoint classes

#include "point.hh"
#include "math.h"

#define MAX_PATH 250

//#define M_PI       3.14159265358979323846
/*
#ifndef sqr
#define sqr(x)                   ((x) * (x))                               // a simple square macro
#endif
*/

#define NUMBER_OF_WORK_STRINGS  20
char    szWork[NUMBER_OF_WORK_STRINGS][MAX_PATH];
int     iWorkPtr = 0;

// Formats the point value for printing, using a static char array
char* iPoint::Format(void)
{
    iWorkPtr++;
    if (iWorkPtr >= NUMBER_OF_WORK_STRINGS)
    {
        iWorkPtr = 0;
    }
    return Format(szWork[iWorkPtr]);
}

// Format using a static string
char* dPoint::Format(void)
{
    iWorkPtr++;
    if (iWorkPtr >= NUMBER_OF_WORK_STRINGS)
    {
        iWorkPtr = 0;
    }
    return Format(szWork[iWorkPtr]);
}

// Format using a static string
char* fPoint::Format(void)
{
    iWorkPtr++;
    if (iWorkPtr >= NUMBER_OF_WORK_STRINGS)
    {
        iWorkPtr = 0;
    }
    return Format(szWork[iWorkPtr]);
}

// returns true if (x,y,x) is between the two values (inclusive)
bool dPoint::IsBetween(dPoint& dpLower, dPoint& dpUpper)
{
    if (x<dpLower.x) return false;
    if (y<dpLower.y) return false;
    if (z<dpLower.z) return false;
    if (x>dpUpper.x) return false;
    if (y>dpUpper.y) return false;
    if (z>dpUpper.z) return false;
    return true;
}

// returns true if (x,y,z) is beteen upper and lower (inclusive)
bool iPoint::IsBetween(iPoint& ipLower, iPoint& ipUpper)
{
    if (x<ipLower.x) return false;
    if (y<ipLower.y) return false;
    if (z<ipLower.z) return false;
    if (x>ipUpper.x) return false;
    if (y>ipUpper.y) return false;
    if (z>ipUpper.z) return false;
    return true;
}

// Calculates the R,T,XYM displacements and stores them in x,y,z respectively
void fPoint::Calculate_R_T_XYM_Disp(fPoint& fpRefPosn, fPoint& fpXYZ_Disp)
{
    // note that fPoint.x is aliased as R
    //           fPoint.y is aliased as Theta
    //           fPoint.z is aliased as X-Y magnitude

    // calculate the new x,y,z position after the displacement
    fPoint  fpNewPosn = fpRefPosn;
    fpNewPosn.Add(fpXYZ_Disp);
    // calculate delta R
    x = (float)(sqrt(fpNewPosn.x * fpNewPosn.x + fpNewPosn.y * fpNewPosn.y) -
                sqrt(fpRefPosn.x * fpRefPosn.x + fpRefPosn.y * fpRefPosn.y) );
    // calculate delta Theta
    double  dRefTheta = atan2(fpRefPosn.y,fpRefPosn.x);
    double  dNewTheta = atan2(fpNewPosn.y,fpNewPosn.x);
    double  dDeltaTheta = dNewTheta - dRefTheta;
    if (dDeltaTheta > M_PI)
    {
        dDeltaTheta -= 2.0 * M_PI;
    }
    else if (dDeltaTheta < -M_PI)
    {
        dDeltaTheta += 2.0 * M_PI;
    }
    y = (float)dDeltaTheta;
    // calculate X-Y magnitude
    double  dDeltaX = fpNewPosn.x - fpRefPosn.x;
    double  dDeltaY = fpNewPosn.y - fpRefPosn.y;
    z = (float)sqrt(dDeltaX * dDeltaX + dDeltaY * dDeltaY);
}

dLine::dLine(void)
{
    dpBase = dPoint(0,0,0);
    dpUnitVector = dPoint(1.0/sqrt(3.0),1.0/sqrt(3.0),1.0/sqrt(3.0));
}

dLine::~dLine(void)
{
}

// construct a dLine from two points
dLine::dLine(dPoint& dpPoint_1, dPoint& dpPoint_2)
{
    dpBase       = dpPoint_1;                       // base point
    dpUnitVector = dpPoint_2;
    dpUnitVector.Subtract(dpPoint_1);               // vector to other point
    dpUnitVector.Divide(dpUnitVector.Magnitude());  // normalize the vector to unit length
}

// returns shortest (i.e., perpendicular) distance from the line to the point
double dLine::DistanceTo(dPoint& dpPoint)
{
    dPoint dpVector(dpPoint);
    dpVector.Subtract(dpBase);      // vector from base point to subject point
    double dDistAlongLine     = dpUnitVector.DotProduct(dpVector);  
    double dDistSqrdBase2Pt   = dpVector.DotProduct(dpVector);
    double dDistLine2Pt       = sqrt(dDistSqrdBase2Pt - sqr(dDistAlongLine));
    return dDistLine2Pt;
}

// considering this dPoint as a vector, compute the dot product with dpVector2
double dPoint::DotProduct(dPoint& dpVector2)
{
    return (this->x*dpVector2.x + 
            this->y*dpVector2.y + 
            this->z*dpVector2.z );
}

// Format a report string for the line
char* dLine::Format(void)
{
    iWorkPtr++;
    if (iWorkPtr >= NUMBER_OF_WORK_STRINGS)
    {
        iWorkPtr = 0;
    }
    sprintf(szWork[iWorkPtr],
            "%s,%s",
            dpBase.Format(),
            dpUnitVector.Format());
    return szWork[iWorkPtr];
}
