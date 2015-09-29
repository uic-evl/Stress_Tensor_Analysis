/* Taken from MTW project by Abedul Haque, December 2009		*/
/* This file has not been modified at all						*/
/* Used for CS3610 class project								*/

#include "Transform.hh"
#include "memory.h"
#include "stdio.h"
#include "assert.h"

// Right multiply this matrix by the matrix tRight and leave the result in tResult
void Transform::RightMult(Transform& tRight, Transform& tResult)
{
    int j,k,l;
    for (j=0; j<4; j++) 
        for (k=0; k<4; k++) 
    {
        tResult.fMatrix[j][k] = 0.0;
        for (l=0; l<4; l++)
        {
            tResult.fMatrix[j][k] += this->fMatrix[j][l] * tRight.fMatrix[l][k];
        }
    }
}

// Right multiply this matrix by the matrix tRight and leave the result in this matrix
void Transform::RightMult(Transform& tRight)
{
    // calculate the matrix product in a temporary area
    Transform tTemp;
    this->RightMult(tRight,tTemp);
    // now copy the matrix product into this object
    int j,k;
    for (j=0; j<4; j++) 
        for (k=0; k<4; k++) 
    {
        this->fMatrix[j][k] = tTemp.fMatrix[j][k];
    }
}

// Left multiply this matrix by the matrix tLeft and leave the result in tResult
void Transform::LeftMult(Transform& tLeft, Transform& tResult)
{
    int j,k,l;
    for (j=0; j<4; j++) 
        for (k=0; k<4; k++) 
    {
        tResult.fMatrix[j][k] = 0.0;
        for (l=0; l<4; l++)
        {
            tResult.fMatrix[j][k] += tLeft.fMatrix[j][l] * this->fMatrix[l][k];
        }
    }
}

// Left multiply this matrix by the matrix tLeft and leave the result in this matrix
//this = tLeft*this
void Transform::LeftMult(Transform& tLeft)
{
    // calculate the matrix product in a temporary area
	
    Transform tTemp;
    this->LeftMult(tLeft,tTemp);
    // now copy the matrix product into this object
    int j,k;
    for (j=0; j<4; j++) 
        for (k=0; k<4; k++) 
    {
        this->fMatrix[j][k] = tTemp.fMatrix[j][k];
    }
}

// Right multiply this matrix by the vector vRight and leave the result in vResult
void Transform::RightMult(Vector& vRight, Vector& vResult)
{
    int j,l;
    for (j=0; j<4; j++) 
    {
        vResult.fVector[j] = 0.0;
        for (l=0; l<4; l++)
        {
            vResult.fVector[j] += this->fMatrix[j][l] * vRight.fVector[l];
        }
    }
    //assert(vResult.VectorW == 1.0f);
}

// Right multiply this matrix by the Vector4 vRight and leave the result in Vector4 vResult
void Transform::RightMult(Vector4& vRight, Vector4& vResult)
{
    int j,l;
    for (j=0; j<4; j++) 
    {
        vResult.fVector[j] = _mm_setzero_ps();
        for (l=0; l<4; l++)
        {
            vResult.fVector[j] = _mm_add_ps(vResult.fVector[j],
                                            _mm_mul_ps(vRight.fVector[l],
                                                       _mm_set_ps1(this->fMatrix[j][l])));
            // vResult.fVector[j] += this->fMatrix[j][l] * vRight.fVector[l];
        }
    }
    assert(1.0f == ((float*)&vResult.fVector[3])[0]);
    assert(1.0f == ((float*)&vResult.fVector[3])[1]);
    assert(1.0f == ((float*)&vResult.fVector[3])[2]);
    assert(1.0f == ((float*)&vResult.fVector[3])[3]);
}

// Left multiply this matrix by the vector vLeft and leave the result in vResult
void Transform::LeftMult(Vector& vLeft, Vector& vResult)
{
    int k,l;
    for (k=0; k<4; k++) 
    {
        vResult.fVector[k] = 0.0;
        for (l=0; l<4; l++)
        {
            vResult.fVector[k] += vLeft.fVector[l] * this->fMatrix[l][k];
        }
    }

    // assert(vResult.VectorW == 1.0f);

    /*
    if (vResult.VectorW == 1.0f) 
    {
        return;
    }   // if (vResult.VectorW == 1.0f) 
    else if (vResult.VectorW != 0)
    {
        for (k=0; k<4; k++) 
        {
            vResult.fVector[k] /= vResult.VectorW;
        }   // for (k=0; k<4; k++) 
        return;
    }   // else if (vResult.VectorW != 0)
    else
    {
        assert(vResult.VectorW == 1.0f);
    }   // else
    */
}

// Left multiply this matrix by the Vector4 vLeft and leave the result in vResult
void Transform::LeftMult(Vector4& vLeft, Vector4& vResult)
{
    int k,l;
    for (k=0; k<4; k++) 
    {
        vResult.fVector[k] = _mm_setzero_ps();
        for (l=0; l<4; l++)
        {
            vResult.fVector[k] = _mm_add_ps(vResult.fVector[k],
                                            _mm_mul_ps(vLeft.fVector[l],
                                                       _mm_set_ps1(this->fMatrix[l][k])));
            // vResult.fVector[k] += vLeft.fVector[l] * this->fMatrix[l][k];
        }
    }
    assert(1.0f == ((float*)&vResult.fVector[3])[0]);
    assert(1.0f == ((float*)&vResult.fVector[3])[1]);
    assert(1.0f == ((float*)&vResult.fVector[3])[2]);
    assert(1.0f == ((float*)&vResult.fVector[3])[3]);
}


// Set up a transformation matrix to translate an object by (ftx,fty,ftz)
void Transform::Translate(float ftx, float fty, float ftz)
{
    memset(fMatrix, 0, 4 * 4 * sizeof(float));
    fMatrix[0][0] = 1.0;
    fMatrix[1][1] = 1.0;
    fMatrix[2][2] = 1.0;
    fMatrix[3][3] = 1.0;
    fMatrix[0][3] = ftx;
    fMatrix[1][3] = fty;
    fMatrix[2][3] = ftz;
}

// Set up a transformation matrix to scale an object by scale factors (fsx,fsy,fsz)
void Transform::Scale(float fsx, float fsy, float fsz)
{
    memset(fMatrix, 0, 4 * 4 * sizeof(float));
    fMatrix[0][0] = fsx;
    fMatrix[1][1] = fsy;
    fMatrix[2][2] = fsz;
    fMatrix[3][3] = 1.0;
}

// Set up a transformation matrix to rotate an object around the X axis by fThetaX radians
void Transform::RotateX(float fThetaX)
{
    memset(fMatrix, 0, 4 * 4 * sizeof(float));
    fMatrix[0][0] = 1.0;
    fMatrix[1][1] = cosf(fThetaX);
    fMatrix[1][2] = -sinf(fThetaX);
    fMatrix[2][2] = cosf(fThetaX);
    fMatrix[2][1] = sinf(fThetaX);
    fMatrix[3][3] = 1.0;
}

// Set up a transformation matrix to rotate an object around the Y axis by fThetaY radians
void Transform::RotateY(float fThetaY)
{
    memset(fMatrix, 0, 4 * 4 * sizeof(float));
    fMatrix[0][0] = cosf(fThetaY);
    fMatrix[0][2] = sinf(fThetaY);
    fMatrix[1][1] = 1.0;
    fMatrix[2][0] = -sinf(fThetaY);
    fMatrix[2][2] = cosf(fThetaY);
    fMatrix[3][3] = 1.0;
}

// Set up a transformation matrix to rotate an object around the Z axis by fThetaZ radians
void Transform::RotateZ(float fThetaZ)
{
    memset(fMatrix, 0, 4 * 4 * sizeof(float));
    fMatrix[0][0] = cosf(fThetaZ);
    fMatrix[0][1] = -sinf(fThetaZ);
    fMatrix[1][0] = sinf(fThetaZ);
    fMatrix[1][1] = cosf(fThetaZ);
    fMatrix[2][2] = 1.0;
    fMatrix[3][3] = 1.0;
}

// Set this matrix equal to the inverse of tInput using Gauss-Jordan elimination with row & column substitution
// This routine was adapted from gaussj() in "Numerical Recipes in C"
void Transform::Invert(Transform& tInput)
{
    int     iColIndex[5],
            iRowIndex[5],
            iPivot   [5],
            i, j, k, l, ll,
            iColumn, iRow;
    double  dWork[5][5],
            dBiggest,
            dDummy,
            dPivotInverse,
            dTemp;

    memset(iColIndex,   0, sizeof(iColIndex));
    memset(iRowIndex,   0, sizeof(iRowIndex));
    memset(iPivot,      0, sizeof(iPivot));
    memset(dWork,       0, sizeof(dWork));

    // Load the source matrix into the double array dWork[][].
    // This is necessary because "Numerical Recipes in C" uses FORTRAN style indexing.
    for (j=1; j<=4; j++) 
        for (i=1; i<=4; i++) 
    {
        dWork[j][i] = tInput.fMatrix[j-1][i-1];
    }   // for (j=1; j<=4; j++) for (i=1; i<=4; i++) 

    // main loop over the columns to be reduced
    for (i=1; i<=4; i++) 
    {
        dBiggest = iRow = iColumn = 0;

        // This is the search for a pivot element
        for (j=1; j<=4; j++) 
            if (iPivot[j] != 1) 
                for (k=1; k<=4; k++) 
        {
            if (iPivot[k] == 0) 
            {
                if (fabs(dWork[j][k]) >= dBiggest) 
                {
                    dBiggest    = fabs(dWork[j][k]);
                    iRow        = j;
                    iColumn     = k;
                }   // if (fabs(dWork[j][k]) >= dBiggest)
            }   // if (iPivot[k] == 0)
        }   // for (j=1; j<=4; j++) if (iPivot[j] != 1) for (k=1; k<=4; k++) 
        assert(dBiggest > 0);

        ++(iPivot[iColumn]);

        // We now have the pivot element, so we interchange rows, if needed, to put the pivot
        // element on the diagonal. The columns are not physically interchanged, only relabeled:
        // iColIndex[i], the column of the ith pivot element, is the ith column that is reduced, while
        // iRowIndex[i] is the row in which that pivot element was originally located. If iRowIndex[i] !=
        // iColIndex[i] there is an implied column interchange. 

        if (iRow != iColumn) 
        {
            for (l=1; l<=4; l++)
            {
                dTemp               = dWork[iRow][l];
                dWork[iRow][l]      = dWork[iColumn][l];
                dWork[iColumn][l]   = dTemp;
            }   // for (l=1; l<=4; l++)
        }   // if (iRow != iColumn) 

        iRowIndex[i] = iRow;    // We are now ready to divide the pivot row by the
        iColIndex[i] = iColumn; // pivot element, located at iRow and iColumn. 

        if (dWork[iColumn][iColumn] == 0.0) 
        {
            printf("\n***** Attempting to invert a singular matrix *****\n");
            QUIT;
        }   // if (dWork[iColumn][iColumn] == 0.0) 

        dPivotInverse           = 1.0 / dWork[iColumn][iColumn];
        dWork[iColumn][iColumn] = 1.0;

        for (l=1;l<=4;l++) 
        {
            dWork[iColumn][l] *= dPivotInverse;
        }   // for (l=1;l<=4;l++) 

        // Next, we reduce the rows... except for the pivot one, of course.
        for (ll=1; ll<=4; ll++) 
            if (ll != iColumn) 
        { 
            dDummy              = dWork[ll][iColumn];
            dWork[ll][iColumn]  = 0.0;
            for (l=1; l<=4; l++) 
            {
                dWork[ll][l]   -= dWork[iColumn][l]*dDummy;
            }   // for (l=1; l<=4; l++) 
        }   // for (ll=1; ll<=4; ll++) if (ll != iColumn) 
    }   // for (i=1; i<=4; i++) 

    // This is the end of the main loop over columns of the reduction. It only remains to unscramble
    // the solution in view of the column interchanges. We do this by interchanging pairs of
    // columns in the reverse order that the permutation was built up.

    for (l=4; l>=1; l--) 
    {
        if (iRowIndex[l] != iColIndex[l]) 
            for (k=1; k<=4; k++)
        {
            dTemp                   = dWork[k][iRowIndex[l]];
            dWork[k][iRowIndex[l]]  = dWork[k][iColIndex[l]];
            dWork[k][iColIndex[l]]  = dTemp;
        }   // if (iRowIndex[l] != iColIndex[l]) for (k=1; k<=4; k++)
    }   // for (l=4; l>=1; l--) 
    
    // Finally, unload the solution matrix back into THIS object's fMatrix[][] array
    for (j=1; j<=4; j++) 
        for (i=1; i<=4; i++) 
    {
        this->fMatrix[j-1][i-1] = (float) dWork[j][i];
    }   // for (j=1; j<=4; j++) for (i=1; i<=4; i++) 

}   // void Transform::Invert(Transform& tInput)


// Print the matrix on the console
void Transform::Print(void)
{
    for (int j=0; j<4; j++)
    {
        printf("%e\t%e\t%e\t%e\n",
            this->fMatrix[j][0],
            this->fMatrix[j][1],
            this->fMatrix[j][2],
            this->fMatrix[j][3]);
    }
}

// Set this matrix equal to the transpose of the tSource matrix
void Transform::Transpose(Transform& tSource)
{
    for (int j=0; j<4; j++) for (int i=0; i<4; i++) 
    {
        this->fMatrix[j][i] = tSource.fMatrix[i][j];
    }
}

// Create a transform to rotate the object around the Z axis, then Y, then X
void Transform::RotateZYX(Vector& vOrientation)
{
    // calculate the three individual transforms
    Transform tXrot,tYrot,tZrot;
    tXrot.RotateX(vOrientation.X());
    tYrot.RotateY(vOrientation.Y());
    tZrot.RotateZ(vOrientation.Z());
    // now calculate THIS transform
    (*this) = tZrot;
    (*this).LeftMult(tYrot);
    (*this).LeftMult(tXrot);
}

/*
// Solve this tranform matrix for the X, Y & Z angles 
Vector& Transform::SolveForAngles(void)
{
// define the subscripts of the general solution matrix
#define CosX_CosZ   0
#define CosX_SinZ   1
#define SinX_CosZ   2
#define SinX_SinZ   3
#define CosX_Index  0
#define SinX_Index  1
#define CosZ_Index  2
#define SinZ_Index  3

    static Vector vOrientation;     // used to return the result
    double dXrot,dYrot,dZrot;
    double dSinY = this->fMatrix[0][2];
    // first handle the general case
    if (fabs(dSinY) > 1e-20)
    {
        // solve for the values of CosX_CosZ, etc.
        Transform tXZcoefficients;
        tXZcoefficients.fMatrix[CosX_CosZ][CosX_CosZ] = 1;              // diagonal
        tXZcoefficients.fMatrix[CosX_SinZ][CosX_SinZ] = 1;              // diagonal
        tXZcoefficients.fMatrix[SinX_CosZ][SinX_CosZ] = 1;              // diagonal
        tXZcoefficients.fMatrix[SinX_SinZ][SinX_SinZ] = 1;              // diagonal
        tXZcoefficients.fMatrix[CosX_CosZ][SinX_SinZ] = -(float)dSinY;  // coefficient
        tXZcoefficients.fMatrix[CosX_SinZ][SinX_CosZ] = +(float)dSinY;  // coefficient
        tXZcoefficients.fMatrix[SinX_CosZ][CosX_SinZ] = +(float)dSinY;  // coefficient
        tXZcoefficients.fMatrix[SinX_SinZ][CosX_CosZ] = -(float)dSinY;  // coefficient
        Transform tXZcoefficientsInverse;
        tXZcoefficientsInverse.Invert(tXZcoefficients);
        Transform tWork = tXZcoefficientsInverse;
        tXZcoefficientsInverse.Transpose(tWork);
        Vector vRightHandSide;
        vRightHandSide.fVector[CosX_CosZ] = this->fMatrix[1][1];
        vRightHandSide.fVector[CosX_SinZ] = this->fMatrix[1][0];
        vRightHandSide.fVector[SinX_CosZ] = this->fMatrix[2][1];
        vRightHandSide.fVector[SinX_SinZ] = this->fMatrix[2][0];
        Vector vSolution;
        tXZcoefficientsInverse.RightMult(vRightHandSide,vSolution);
        // OK, we have CosX_CosZ, etc.  Now solve for CosX, etc.
        Transform tXZcoefficients_2;
        tXZcoefficients_2.fMatrix[CosX_Index][CosX_Index] = (float)(+1.0 / vSolution.fVector[CosX_CosZ]);    // diagonal
        tXZcoefficients_2.fMatrix[CosX_Index][SinX_Index] = (float)(-1.0 / vSolution.fVector[SinX_CosZ]);    // coefficient
        tXZcoefficients_2.fMatrix[SinX_Index][SinX_Index] = (float)(+1.0 / vSolution.fVector[SinX_SinZ]);    // diagonal
        tXZcoefficients_2.fMatrix[SinX_Index][CosX_Index] = (float)(-1.0 / vSolution.fVector[CosX_SinZ]);    // coefficient
        tXZcoefficients_2.fMatrix[CosZ_Index][CosZ_Index] = (float)(+1.0 / vSolution.fVector[CosX_CosZ]);    // diagonal
        tXZcoefficients_2.fMatrix[CosZ_Index][SinZ_Index] = (float)(-1.0 / vSolution.fVector[CosX_SinZ]);    // coefficient
        tXZcoefficients_2.fMatrix[SinZ_Index][SinZ_Index] = (float)(+1.0 / vSolution.fVector[SinX_SinZ]);    // diagonal
        tXZcoefficients_2.fMatrix[SinZ_Index][CosZ_Index] = (float)(-1.0 / vSolution.fVector[SinX_CosZ]);    // coefficient
        Transform tXZcoefficientsInverse_2;
        tXZcoefficientsInverse_2.Invert(tXZcoefficients_2);
        Vector vRightHandSide_2,vSolution_2;
        memset(&vRightHandSide_2,0,sizeof(vRightHandSide_2));
        tXZcoefficientsInverse_2.RightMult(vRightHandSide_2,vSolution_2);
        // OK, we now have CosX, etc.
        dXrot = (float)atan2(vSolution_2.fVector[SinX_Index],vSolution_2.fVector[CosX_Index]);
        dYrot = (float)atan2(dSinY,  this->fMatrix[0][0] / vSolution_2.fVector[CosZ_Index]);
        dZrot = (float)atan2(vSolution_2.fVector[SinZ_Index],vSolution_2.fVector[CosZ_Index]);
        vOrientation = Vector(dXrot,dYrot,dZrot);
    }
    else // (dSinY == 0)
    {
    }
    return vOrientation;
}
*/

//*
// Solve this transform matrix for the X, Y & Z angles using Scott Tashman's dcm_to_rot()
Vector& Transform::SolveForAngles(void)
{
    double  dDirectionCosineMatrix[3][3];           // used to pass DCM into dcm_to_rot()
    double  dXrotation, dYrotation, dZrotation;     // used to receive results from dcm_to_rot()
    static Vector vOrientation;                     // used to return results to the calling routine

    // copy the DCM out of THIS transform and into dDirectionCosineMatrix[][]
    for (int j=0; j<3; j++) for (int i=0; i<3; i++) 
    {
        dDirectionCosineMatrix[j][i] = this->fMatrix[i][j];
    }   // for (int j=0; j<3; j++) for (int i=0; i<3; i++) 

    // elicit the solution from dcm_to_rot()
    dcm_to_rot( dDirectionCosineMatrix,
                RT_BODY,
                3,
                2,
                1,
                &dZrotation,
                &dYrotation,
                &dXrotation);

    // map all angles to the range of -M_PI => +M_PI
    if (dXrotation < -M_PI) dXrotation += 2.0 * M_PI;
    if (dYrotation < -M_PI) dYrotation += 2.0 * M_PI;
    if (dZrotation < -M_PI) dZrotation += 2.0 * M_PI;
    if (dXrotation > +M_PI) dXrotation -= 2.0 * M_PI;
    if (dYrotation > +M_PI) dYrotation -= 2.0 * M_PI;
    if (dZrotation > +M_PI) dZrotation -= 2.0 * M_PI;

    return vOrientation = new Vector(-dXrotation, -dYrotation, -dZrotation);
}
//*/

/*
// Solve this tranform matrix for the X, Y & Z angles 
Vector& Transform::SolveForAngles(void)
{
    static Vector vOrientation;     // used to return the result
    double dSinY = this->fMatrix[0][2];
    double dCosY = sqrt(-(fMatrix[2][2] * fMatrix[0][2] * fMatrix[0][1] + fMatrix[1][2] * fMatrix[0][0]) / fMatrix[2][1]);
    double dSinX = -fMatrix[1][2] / dCosY;
    double dCosX = +fMatrix[2][2] / dCosY;
    double dSinZ = -fMatrix[0][1] / dCosY;
    double dCosZ = +fMatrix[0][0] / dCosY;
    vOrientation = Vector(atan2(dSinX,dCosX),
                          atan2(dSinY,dCosY),
                          atan2(dSinZ,dCosZ));
    return vOrientation;
}
//*/