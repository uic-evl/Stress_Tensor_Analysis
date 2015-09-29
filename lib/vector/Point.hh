// point.h -- header for point.cpp

#pragma once

#include "stdio.h"
#include "math.h"
#include "assert.h"

#ifndef sqr
#define sqr(x)                   ((x) * (x))                               // a simple square macro
#endif
#ifndef QUIT
#define QUIT {printf("\n***** Fatal error on line %d of '%s'\n",__LINE__,__FILE__);assert(0);exit(999);} // fatal-error macro
#endif
#ifndef NOALLOC
#define NOALLOC(x) {printf("\n***** Unable to allocate enough memory for %s\n",x);QUIT;}
#endif

class iPoint
{
public:
    int x,y,z;  // coordinates of point

    iPoint(void)
    {
        x = 0;
        y = 0;
        z = 0;
    };
    iPoint(int ix, int iy, int iz)
    {
        x = ix;
        y = iy;
        z = iz;
    };
    iPoint(double fx, double fy, double fz)
    {
        x = (int) fx;
        y = (int) fy;
        z = (int) fz;
    };
    void    Add(iPoint& ip)
    {
        x += ip.x;
        y += ip.y;
        z += ip.z;
    };
    void    Subtract(iPoint& ip)
    {
        x -= ip.x;
        y -= ip.y;
        z -= ip.z;
    };
    void    Multiply(int iMultiplier)
    {
        x *= iMultiplier;
        y *= iMultiplier;
        z *= iMultiplier;
    };
    void    Divide(int iDivisor)
    {
        x /= iDivisor;
        y /= iDivisor;
        z /= iDivisor;
    };
    char*   Format(char* szBfr)
    {
        sprintf(szBfr,"(%3d,%3d,%3d)",x,y,z);
        return szBfr;
    };

    bool IsEqualTo(iPoint& ipComp)
    {
        if (ipComp.x == x && ipComp.y == y && ipComp.z == z)
        {
            return true;
        }
        else
        {
            return false;
        }
    };

    // Calculates the RMS magnitude of the vector (rounded)
    double Magnitude(void)
    {
        return (sqrt((double)x*x + (double)y*y + (double)z*z) + 0.5);
    }

    // Formats the point value for printing, using a static char array
    char* Format(void);

    // calculate the product of the three coordinates
    int Volume(void)
    {
        return x*y*z;
    }
    // returns true if (x,y,z) is beteen upper and lower (inclusive)
    bool IsBetween(iPoint& ipLower, iPoint& ipUpper);
};

class fPoint
{
public:
    // X value of vector; R value in polar coordinates
    float x;
    // Y value of vector; Theta value in polar coordinates
    float y;
    // Z value of vector; X-Y magnitude in polar coordinates
    float z;

    fPoint(void)
        : x(0)
    {
        x = 0.0f;
        y = 0.0f;
        z = 0.0f;
    };
    fPoint(float ix, float iy, float iz)
    {
        x = ix;
        y = iy;
        z = iz;
    };
    void    Add(iPoint& ip)
    {
        x += ip.x;
        y += ip.y;
        z += ip.z;
    };
    void    Add(fPoint& ip)
    {
        x += ip.x;
        y += ip.y;
        z += ip.z;
    };
    void    Subtract(iPoint& ip)
    {
        x -= ip.x;
        y -= ip.y;
        z -= ip.z;
    };
    void    Multiply(int iMultiplier)
    {
        x *= iMultiplier;
        y *= iMultiplier;
        z *= iMultiplier;
    };
    void    Divide(int iDivisor)
    {
        x /= iDivisor;
        y /= iDivisor;
        z /= iDivisor;
    };
    void    Subtract(fPoint& ip)
    {
        x -= ip.x;
        y -= ip.y;
        z -= ip.z;
    };
    void    Multiply(float iMultiplier)
    {
        x *= iMultiplier;
        y *= iMultiplier;
        z *= iMultiplier;
    };
    void    Divide(float iDivisor)
    {
        x /= iDivisor;
        y /= iDivisor;
        z /= iDivisor;
    };
    char* Format(char* szBfr)
    {
        sprintf(szBfr,"(%6.3f,%6.3f,%6.3f)",x,y,z);
        return szBfr;
    };

    void Round(void)
    {
        x = (float)floor(x + 0.5f);
        y = (float)floor(y + 0.5f);
        z = (float)floor(z + 0.5f);
    };

    bool IsZero(void)
    {
        if (x || y || z) return false;
        else             return true;
    };

    // Calculates the RMS magnitude of the vector
    float Magnitude(void)
    {
        return (float)sqrt(x*x + y*y + z*z);
    }

    // Format using a static string
    char* Format(void);

    // overloaded constructor for fPoint
    fPoint(iPoint& ipValue)
    {
        x = (float)ipValue.x;
        y = (float)ipValue.y;
        z = (float)ipValue.z;
    }

    // compare this fPoint to another one
    bool IsEqual(fPoint& fpArg)
    {
        if (x != fpArg.x) return false;
        if (y != fpArg.y) return false;
        if (z != fpArg.z) return false;
        return true;
    }
    // Calculates the R,T,XYM displacements and stores them in x,y,z respectively
    void Calculate_R_T_XYM_Disp(fPoint& fpRefPosn, fPoint& fpXYZ_Disp);

    // construct from a dPoint object
    fPoint(double dx, double dy, double dz)
    {
        x = (float) dx;
        y = (float) dy;
        z = (float) dz;
    }
    // construct fPoint object from integer data
    fPoint(int ix, int iy, int iz)
    {
        x = (float) ix;
        y = (float) iy;
        z = (float) iz;
    }
};

class dPoint
{
public:
    double x,y,z;  // coordinates of point

    dPoint(void)
    {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    };
    dPoint(double ix, double iy, double iz)
    {
        x = ix;
        y = iy;
        z = iz;
    };
    void    Add(iPoint& ip)
    {
        x += ip.x;
        y += ip.y;
        z += ip.z;
    };
    void    Add(fPoint& ip)
    {
        x += ip.x;
        y += ip.y;
        z += ip.z;
    };
    void    Add(dPoint& ip)
    {
        x += ip.x;
        y += ip.y;
        z += ip.z;
    };
    void    Subtract(iPoint& ip)
    {
        x -= ip.x;
        y -= ip.y;
        z -= ip.z;
    };
    void    Subtract(fPoint& ip)
    {
        x -= ip.x;
        y -= ip.y;
        z -= ip.z;
    };
    void    Subtract(dPoint& ip)
    {
        x -= ip.x;
        y -= ip.y;
        z -= ip.z;
    };
    void    Multiply(int iMultiplier)
    {
        x *= iMultiplier;
        y *= iMultiplier;
        z *= iMultiplier;
    };
    void    Multiply(double iMultiplier)
    {
        x *= iMultiplier;
        y *= iMultiplier;
        z *= iMultiplier;
    };
    void    Divide(int iDivisor)
    {
        x /= iDivisor;
        y /= iDivisor;
        z /= iDivisor;
    };
    void    Divide(double iDivisor)
    {
        x /= iDivisor;
        y /= iDivisor;
        z /= iDivisor;
    };
    char* Format(char* szBfr)
    {
        sprintf(szBfr,"(%10.6f,%10.6f,%10.6f)",x,y,z);
        return szBfr;
    };
    void Round(void)
    {
        x = floor(x + 0.5);
        y = floor(y + 0.5);
        z = floor(z + 0.5);
    }

    // // normalize this vector to a length of 1.0
    void Normalize(void)
    {
        double dLength = sqrt(x*x + y*y + z*z);
        // beware of normalizing very small vectors
        if (dLength > 1.0e-10)
        {
            Divide(dLength);
        }
        else
        {
            // if the vector is too small to normalize,
            // we will just assign an arbitrary value.
            // we do this so as to not screw up the quadratic search
            x = y = z = (1.0 / sqrt(3.0));
        }
    }

    // Return true if dpTest is identical to this object
    bool IsEqual(dPoint& dpTest)
    {
        if (x == dpTest.x &&
            y == dpTest.y &&
            z == dpTest.z)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    // moves the base point in the direction of dpGradient for the distance of dDistance
    void GradientMove(dPoint& dpGradient, double dDistance)
    {
        dPoint  dpDelta(dpGradient);
        dpDelta.Multiply(dDistance);
        x += dpDelta.x;
        y += dpDelta.y;
        z += dpDelta.z;
    }

    // Calculates the RMS magnitude of the vector
    double Magnitude(void)
    {
        return sqrt(x*x + y*y + z*z);
    }
    // Format using a static string
    char* Format(void);
    // returns true if (x,y,x) is between the two values (inclusive)
    bool IsBetween(dPoint& dpLower, dPoint& dpUpper);
    // considiering this dPoint as a vector, compute the dot product with dpVector2
    double DotProduct(dPoint& dpVector2);
};

class dLine
{
public:
    dLine(void);
    ~dLine(void);
    // a point on the line
    dPoint dpBase;
    // a vector of unit length that establishes the direction of the line
    dPoint dpUnitVector;
    // construct a dLine from two points
    dLine(dPoint& dpPoint_1, dPoint& dpPoint_2);
    // returns shortest (i.e., perpendicular) distance from the line to the point
    double DistanceTo(dPoint& dpPoint);
    // Format a report string for the line
    char* Format(void);
};

class iArea{

public:
	iPoint ipMin ;
	iPoint ipMax ;
	bool bIsSet ;

	iArea(void) { bIsSet = false ; }

	iArea( iArea& iaROI){
		this->ipMin.x = iaROI.ipMin.x ;
		this->ipMin.y = iaROI.ipMin.y ;

		this->ipMax.x = iaROI.ipMax.x ;
		this->ipMax.y = iaROI.ipMax.y ;
		this->bIsSet = true ;
	}

	iArea(iPoint& tl, iPoint& br ){
		this->ipMin.x = tl.x ;
		this->ipMin.y = tl.y ;
		//this->ipMin.z = tl.z ;
		this->ipMax.x = br.x ;
		this->ipMax.y = br.y ;
		//this->ipMax.z = br.z ;
		bIsSet = true ;
	}

	iArea( int x1, int y1, int x2, int y2){
		this->ipMin.x = x1 ;
		this->ipMin.y = y1 ;
		//this->ipMin.z = z1 ;
		this->ipMax.x = x2 ;
		this->ipMax.y = y2 ;
		//this->ipMax.z = z2 ;
		this->bIsSet = true ;
	}
	~iArea(void){} ;
	
} ;