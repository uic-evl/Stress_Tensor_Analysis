// Vector.hh -- Class definition for Vector, a class of homogenous vectors

#pragma once

#define  _USE_MATH_DEFINES

#include "math.h"
#include "stdio.h"
#include "assert.h"
#include "emmintrin.h"      // header file for SSE intrinsics

#define VectorX fVector[0]
#define VectorY fVector[1]
#define VectorZ fVector[2]
#define VectorW fVector[3]

#ifndef sqr
#define sqr(x) ((x)*(x))
#endif

#define MAX_PATH 250

#ifndef min
#define min(a,b) (a<b?a:b)
#endif

#ifndef max
#define max(a,b) (a>b?a:b)
#endif

class Vector
{
public:
    float fVector[4];
    Vector(void){}
    ~Vector(void){}

    // Constructor with three float values
    Vector(float x, float y, float z)
    {
        VectorX = x;
        VectorY = y;
        VectorZ = z;
        VectorW = 1.0;
    }

    // Constructor for three integer values
    Vector(int ix, int iy, int iz)
    {
        VectorX = (float)ix;
        VectorY = (float)iy;
        VectorZ = (float)iz;
        VectorW = 1.0;
    }

    // Constructor for an existing Vector
    Vector(Vector& v)
    {
        VectorX = v.VectorX;
        VectorY = v.VectorY;
        VectorZ = v.VectorZ;
        VectorW = v.VectorW;
    }

    // Constructor for a pointer to an existing Vector
    Vector(Vector* pv)
    {
        VectorX = pv->VectorX;
        VectorY = pv->VectorY;
        VectorZ = pv->VectorZ;
        VectorW = pv->VectorW;
    }

    // Constructor for three double values
    Vector(double dx, double dy, double dz)
    {
        VectorX = (float)dx;
        VectorY = (float)dy;
        VectorZ = (float)dz;
        VectorW = 1.0;
    }

    // Extract X coordinate
    float X(void)
    {
        //assert(VectorW == 1.0f);
        return VectorX;
    }

    // Extract Y coordinate
    float Y(void)
    {
        //assert(VectorW == 1.0f);
        return VectorY;
    }

    // Extract Z coordinate
    float Z(void)
    {
        //assert(VectorW == 1.0f);
        return VectorZ;
    }

    // Calculate the magnitude of the vector
    float Magnitude(void)
    {
        return (float)sqrt(sqr(X()) + 
                           sqr(Y()) + 
                           sqr(Z()) );
    }

    // Calculate dot product with another vector
    float Dot(Vector& vDot)
    {
        return X() * vDot.X() + 
               Y() * vDot.Y() + 
               Z() * vDot.Z() ;
    }

    // Dot product with pointer to other vector
    float Dot(Vector* vp)
    {
        return Dot(*vp);
    }

    // Normalize the vector to a length of 1.0
    void Normalize(void)
    {
        float fMag = Magnitude();
        if (fMag != 0)
        {
            VectorX /= fMag;
            VectorY /= fMag;
            VectorZ /= fMag;
            VectorW  = 1;
        }
        else
        {
            VectorX = 0;
            VectorY = 0;
            VectorZ = 0;
            VectorW = 1;
        }
    }

    // Calculate cross product with b and leave in vCross
    void Cross(Vector& vb, Vector& vCross)
    {
        vCross.VectorX = Y() * vb.Z() - Z() * vb.Y();
        vCross.VectorY = Z() * vb.X() - X() * vb.Z();
        vCross.VectorZ = X() * vb.Y() - Y() * vb.X();
        vCross.VectorW = 1.0;
    }

    // Compute cross product using pointers
    void Cross(Vector* pb, Vector* pCross)
    {
        Cross(*pb,*pCross);
    }

    // Multiply the vector by the scalar value fMult
    void ScalarMult(float fMult)
    {
        VectorX *= fMult;
        VectorY *= fMult;
        VectorZ *= fMult;
    }

    // Multiply this vector by -1.0
    void Negate(void)
    {
        ScalarMult(-1.0);
    }

    // Print the vector on the console
    void Print(void)
    {
        printf("%10.3f\t%10.3f\t%10.3f\n",
            this->X(),
            this->Y(),
            this->Z());
    }

    // Set each coordinate of this vector to the minimum of this vector and vSource
    void Min(Vector& vSource)
    {
        this->VectorX = min(this->X(),vSource.X());
        this->VectorY = min(this->Y(),vSource.Y());
        this->VectorZ = min(this->Z(),vSource.Z());
    }

    // Set each coordinate of this vector to the maximum of this vector and vSource
    void Max(Vector& vSource)
    {
        this->VectorX = max(this->X(),vSource.X());
        this->VectorY = max(this->Y(),vSource.Y());
        this->VectorZ = max(this->Z(),vSource.Z());
    }

    // Subtract the vSubtrahend coordinates from this object
    void Subtract(Vector& vSubtrahend)
    {
        this->VectorX = this->X() - vSubtrahend.X();
        this->VectorY = this->Y() - vSubtrahend.Y();
        this->VectorZ = this->Z() - vSubtrahend.Z();
    }

    // Calculate the floor() function on each coordinate of this vector
    void Floor(void)
    {
        this->VectorX = (float)floor(this->X());
        this->VectorY = (float)floor(this->Y());
        this->VectorZ = (float)floor(this->Z());
    }

    // Calculate the ceil() function on each coordinate of this vector
    void Ceil(void)
    {
        this->VectorX = (float)ceil(this->X());
        this->VectorY = (float)ceil(this->Y());
        this->VectorZ = (float)ceil(this->Z());
    }

    // Returns a pointer to a string with the vector in ASCII.
    char* Format(void)
    {
        static char szResult[MAX_PATH];
        sprintf(szResult,"(%6.3f,%6.3f,%6.3f)",X(),Y(),Z());
        return szResult;
    }

    // Add the Vector to this one
    void Add(Vector& vAddend)
    {
        this->VectorX = this->X() + vAddend.X();
        this->VectorY = this->Y() + vAddend.Y();
        this->VectorZ = this->Z() + vAddend.Z();
    }

    // Convert this Vector from degrees to radians
    Vector& Degrees_to_Radians(void);
    // Convert this Vector from radians to degrees
    Vector& Radians_to_Degrees(void);

    // Static function returns the null vector (0,0,0).
    static Vector& Null(void)
    {
        // format a null vector
        static Vector vNull(0,0,0);
        return vNull;
    }

    // Calculates ((THIS cross Yaxis) dot Zaxis); positive result means right-hand system, negative result means left-hand.
    float CrossDot(Vector& Yaxis, Vector& Zaxis)
    {
        Vector vCross;
        this->Cross(Yaxis, vCross);
        float fCrossDot = vCross.Dot(Zaxis);
        return fCrossDot;
    }

	// Interpolation of two vectors. fDelta should be between 0 to 1
	// if fDelta is 0, the resulting vector is current vector. if fDelta is 1, the resulting vector is vEnd
	void Interpolate(Vector& vEnd, float fDelta){
		this->VectorX = (1.0 - fDelta)*this->VectorX + fDelta*vEnd.VectorX ; 
		this->VectorY = (1.0 - fDelta)*this->VectorY + fDelta*vEnd.VectorY ; 
		this->VectorZ = (1.0 - fDelta)*this->VectorZ + fDelta*vEnd.VectorZ ; 
        this->VectorW = 0;		
		//this->Normalize() ; 
	}

};

// this class packages four values in each
class Vector4
{
public:
    __m128 fVector[4];
    Vector4(void){}
    ~Vector4(void){}

    // Constructor with three float values
    Vector4(float x, float y, float z)
    {
        VectorX = _mm_set_ps1(x);
        VectorY = _mm_set_ps1(y);
        VectorZ = _mm_set_ps1(z);
        VectorW = _mm_set_ps1(1.0f);
    }
    // constructor from another Vector4
    Vector4(const Vector4& v4Source)
    {
        this->fVector[0] = v4Source.fVector[0];
        this->fVector[1] = v4Source.fVector[1];
        this->fVector[2] = v4Source.fVector[2];
        this->fVector[3] = v4Source.fVector[3];
    }
	Vector4(Vector4* v4Source)
    {
        this->fVector[0] = v4Source->fVector[0];
        this->fVector[1] = v4Source->fVector[1];
        this->fVector[2] = v4Source->fVector[2];
        this->fVector[3] = v4Source->fVector[3];
    }

	void Assign(float x, float y, float z){
		
		VectorX = _mm_set_ps1(x);
        VectorY = _mm_set_ps1(y);
        VectorZ = _mm_set_ps1(z);
        VectorW = _mm_set_ps1(1.0f);	
	}
    // Multiply the vector by the scalar value fMult
    void ScalarMult(float fMult)
    {
        static __m128 m128_fMult;
        m128_fMult = _mm_set_ps1(fMult);
        VectorX = _mm_mul_ps(VectorX,m128_fMult);
        VectorY = _mm_mul_ps(VectorY,m128_fMult);
        VectorZ = _mm_mul_ps(VectorZ,m128_fMult);
    }
    // Multiply the vector by the four scalar values in Vector fMult
    void ScalarMult(Vector& fMult)
    {
        __m128* m128_fMult = (__m128*)&fMult;       // recast the "Vector" as an __m128
        VectorX = _mm_mul_ps(VectorX,*m128_fMult);
        VectorY = _mm_mul_ps(VectorY,*m128_fMult);
        VectorZ = _mm_mul_ps(VectorZ,*m128_fMult);
    }
    // Multiply the four vectors by the four scalar values in m128_fMult
    void ScalarMult(__m128* m128_fMult)
    {
        this->fVector[0] = _mm_mul_ps(this->fVector[0],*m128_fMult);
        this->fVector[1] = _mm_mul_ps(this->fVector[1],*m128_fMult);
        this->fVector[2] = _mm_mul_ps(this->fVector[2],*m128_fMult);
    }
    // Add the Vector4 vAddend to this object
    void Add(Vector4& vAddend)
    {
        this->VectorX = _mm_add_ps(this->VectorX,vAddend.VectorX);
        this->VectorY = _mm_add_ps(this->VectorY,vAddend.VectorY);
        this->VectorZ = _mm_add_ps(this->VectorZ,vAddend.VectorZ);
        this->VectorW = _mm_set_ps1(1.0f);
    }
    // Add the same Vector vAddend to each vector of this Vector4 object
    void Add(Vector& vAddend)
    {
        this->fVector[0] = _mm_add_ps(this->fVector[0],_mm_set_ps1(vAddend.fVector[0]));
        this->fVector[1] = _mm_add_ps(this->fVector[1],_mm_set_ps1(vAddend.fVector[1]));
        this->fVector[2] = _mm_add_ps(this->fVector[2],_mm_set_ps1(vAddend.fVector[2]));
        this->fVector[3] = _mm_set_ps1(1.0f);
    }
    // Subtract the Vector4 vSubtrahend from this object
    void Subtract(Vector4& vSubtrahend)
    {
        this->VectorX = _mm_sub_ps(this->VectorX,vSubtrahend.VectorX);
        this->VectorY = _mm_sub_ps(this->VectorY,vSubtrahend.VectorY);
        this->VectorZ = _mm_sub_ps(this->VectorZ,vSubtrahend.VectorZ);
        this->VectorW = _mm_set_ps1(1.0f);
    }
    // Subtract the same Vector vSubtrahend from each of the four vectors in this object
    void Subtract(Vector& vSubtrahend)
    {
        this->VectorX = _mm_sub_ps(this->VectorX,_mm_set_ps1(vSubtrahend.VectorX));
        this->VectorY = _mm_sub_ps(this->VectorY,_mm_set_ps1(vSubtrahend.VectorY));
        this->VectorZ = _mm_sub_ps(this->VectorZ,_mm_set_ps1(vSubtrahend.VectorZ));
        this->VectorW = _mm_set_ps1(1.0f);
    }
    // Returns a pointer to an __m128 with the magnitudes of the four vectors in this Vector4
    Vector* Magnitude(void)
    {
        static __m128 m128_Result;
        m128_Result = _mm_mul_ps(this->VectorX,this->VectorX);
        m128_Result = _mm_add_ps(m128_Result,_mm_mul_ps(this->VectorY,this->VectorY));
        m128_Result = _mm_add_ps(m128_Result,_mm_mul_ps(this->VectorZ,this->VectorZ));
        m128_Result = _mm_sqrt_ps(m128_Result);
        return (Vector*) &m128_Result;
    }
    // Returns a pointer to an __m128 with the four distances between the four points in this Vector4 and those in the argument
    Vector* Distance(Vector4& v4toPoints)
    {
        static __m128 m128_Result;

        // calculate the differences by component
        __m128 m128_VectorX = _mm_sub_ps(this->VectorX, v4toPoints.VectorX);
        __m128 m128_VectorY = _mm_sub_ps(this->VectorY, v4toPoints.VectorY);
        __m128 m128_VectorZ = _mm_sub_ps(this->VectorZ, v4toPoints.VectorZ);

        // accumulate the square of the differences
        m128_Result         = _mm_add_ps(_mm_add_ps(_mm_mul_ps(m128_VectorX,  m128_VectorX),
                                                    _mm_mul_ps(m128_VectorY,  m128_VectorY)),
                                                    _mm_mul_ps(m128_VectorZ,  m128_VectorZ) );

        // take the square root of the accumulations
        m128_Result         = _mm_sqrt_ps(m128_Result);

        return (Vector*) &m128_Result;
    }
};
