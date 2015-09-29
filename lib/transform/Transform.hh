/* Taken from MTW project by Abedul Haque, December 2009		*/
/* This file has not been modified at all						*/
/* used for CS3610 class project								*/

#pragma once

#include "../vector/Vector.hh"
#include "../math/MathTools.hh"
#include <string.h>

// standard macros
#ifndef sqr
#define sqr(x)                   ((x) * (x))                               // a simple square macro
#endif
#ifndef QUIT
#define QUIT {printf("\n***** Fatal error on line %d of '%s'\n",__LINE__,__FILE__);assert(0);exit(999);} // fatal-error macro
#endif
#ifndef NOALLOC
#define NOALLOC(x) {printf("\n***** Unable to allocate enough memory for %s\n",x);QUIT;}
#endif


class Transform
{
public: // member variables
    float fMatrix[4][4];

public: // member functions
    Transform(void)
    {
        memset(this,0,(4 * 4 * sizeof(float)));
    }
    ~Transform(void){}
    // Right multiply this matrix by the matrix tRight and leave the result in tResult
    void RightMult(Transform& tRight, Transform& tResult);
    // Right multiply this matrix by the matrix *tRight and leave the result in *tResult (pointer version)
    void RightMult(Transform *tRight, Transform *tResult)
    {
        RightMult(*tRight, *tResult);
    }
    // Right multiply this matrix by the matrix tRight and leave the result in this matrix
    void RightMult(Transform& tRight);
    // Right multiply this matrix by the matrix tRight and leave the result in this matrix (pointer version)
    void RightMult(Transform *tRight)
    {
        RightMult(*tRight);
    }
    // Left multiply this matrix by the matrix tLeft and leave the result in tResult
    void LeftMult(Transform& tLeft, Transform& tResult);
    // Left multiply this matrix by the matrix *tLeft and leave the result in *tResult (pointer version)
    void LeftMult(Transform *tLeft, Transform *tResult)
    {
        LeftMult(*tLeft, *tResult);
    }
    // Left multiply this matrix by the matrix tLeft and leave the result in this matrix
    void LeftMult(Transform& tLeft);
    // Left multiply this matrix by the matrix tLeft and leave the result in this matrix (pointer version)
    void LeftMult(Transform *tLeft)
    {
        LeftMult(*tLeft);
    }
    // Right multiply this matrix by the vector vRight and leave the result in vResult
    void RightMult(Vector& vRight, Vector& vResult);
    // Right multiply this matrix by the vector vRight and leave the result in vResult (pointer version)
    void RightMult(Vector *pRight, Vector *pResult)
    {
        RightMult(*pRight, *pResult);
    }
    // Right multiply this matrix by the Vecto4 vRight and leave the result in Vector4 vResult
    void RightMult(Vector4& vRight, Vector4& vResult);
    // Left multiply this matrix by the vector vLeft and leave the result in vResult
    void LeftMult(Vector& vLeft, Vector& vResult);
    // Left multiply this matrix by the Vector4 vLeft and leave the result in vResult
    void LeftMult(Vector4& vLeft, Vector4& vResult);
    // Left multiply this matrix by the vector vLeft and leave the result in vResult (pointer version)
    void LeftMult(Vector *pLeft, Vector *pResult)
    {
        LeftMult(*pLeft, *pResult);
    }
    // Set up a transformation matrix to translate the point by (ftx,fty,ftz)
    void Translate(float ftx, float fty, float ftz);
    // Set up a transformation matrix to translate the point by Vector v
    void Translate(Vector& v)
    {
        Translate(v.X(), v.Y(), v.Z());
    }
    // Set up a transformation matrix to translate from point v to origin
    void MoveToOrigin(Vector& v)
    {
        Translate(-v.X(), -v.Y(), -v.Z());
    }
    // Set up a transformation matrix to scale by factors (fsx,fsy,fsz)
    void Scale(float fsx, float fsy, float fsz);
    // Set up a transformation matrix to scale with the coordinates of vScale
    void Scale(Vector& vScale)
    {
        Scale(vScale.X(), vScale.Y(), vScale.Z());
    }
    // Set up a transformation matrix to rotate the object around the X axis by fThetaX radians
    void RotateX(float fThetaX);
    // Set up a transformation matrix to rotate an object around the Y axis by fThetaY radians
    void RotateY(float fThetaY);
    // Set up a transformation matrix to rotate an object around the Z axis by fThetaZ radians
    void RotateZ(float fThetaZ);
    // Set up a transformation matrix to rotate the object around the X axis by dThetaX radians
    void RotateX(double dThetaX)
    {
        RotateX((float)dThetaX);
    }
    // Set up a transformation matrix to rotate an object around the Y axis by dThetaY radians
    void RotateY(double dThetaY)
    {
        RotateY((float)dThetaY);
    }
    // Set up a transformation matrix to rotate an object around the Z axis by dThetaZ radians
    void RotateZ(double dThetaZ)
    {
        RotateZ((float)dThetaZ);
    }
    // Set this matrix equal to the inverse of tInput using Gauss-Jordan elimination with row & column substitution
    void Invert(Transform& tInput);
    // Print the matrix on the console
    void Print(void);
    // Set this matrix equal to the transpose of the tSource matrix
    void Transpose(Transform& tSource);
    // Create a transform to rotate the object around the Z axis, then Y, then X
    void RotateZYX(Vector& vOrientation);
    // Solve THIS transform matrix for the X, Y & Z angles using Scott Tashman's dcm_to_rot()
    Vector& SolveForAngles(void);

    // Static function returns the identity traansform matrix
    static Transform& Identity(void)
    {
        // format an identity transformation matrix
        static Transform tIdentity;
        memset(&tIdentity,0,sizeof(tIdentity));
        tIdentity.fMatrix[0][0] = 1;
        tIdentity.fMatrix[1][1] = 1;
        tIdentity.fMatrix[2][2] = 1;
        tIdentity.fMatrix[3][3] = 1;
        return tIdentity;
    }

};
