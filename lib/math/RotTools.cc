/*
  rottools - created 3/11/91 (S. Tashman)

  routines for dealing with rigid body rotations in 3-D space.

  dcm_to_rot  - determines ordered rotation angles from dir cos matrix
  rot_to_dcm  - determines dir cos matrix from ordered rotation angles
  orthonorm   - creates an orthonormal basis from three noncoplanar vectors

  These routines also use routines defined in "mathtools.c".
*/

// #ifndef SuppressStandardHeaders
//#include "stdafx.h"
// #endif

#include <math.h>
#include <stdio.h>
#include "MathTools.hh"

/***********************************************************************
  dcm_to_rot - takes a direction cos matrix a[3][3] and calculates
  body- or space-fixed rotation angles r1,r2,r3 in the order defined 
  by d1,d2,d3

  NOTE: this algorithm always assumes that  -pi/2 < r2 < pi/2, and
  calculates a solution accordingly.

  This algorithm only good for 3 different/independent rotations 
  (i.e. 1,2,3), not when performing rotations such as 2,1,2.

*/

/*
#define RT_FWD      1
#define RT_REV     -1
#define RT_SPACE    1
#define RT_BODY    -1
*/

void dcm_to_rot(double  a[3][3], 
                int     type, 
                int     d1, 
                int     d2, 
                int     d3, 
                double* r1, 
                double* r2, 
                double* r3)

// double a[3][3];			/* input dir cos matrix */
// int type;			    /* 1 for space-fixed, -1 for body-fixed */
// int d1, d2, d3;			/* order for rotations */
// double *r1, *r2, *r3;    /* pointers to returned rotations */

{
  double b[3][3];   /* temp storage for rearranged el's of a */
  int d[3];			/* temp storage for --d1,--d2,--d3 */
  double sign;
  double alpha;
  double beta;
  double deter;
  double PI;
  int i, j;
  int bi, bj;
  int rotdir;

  PI = 4.0*atan(1.0);

/* check validity of input matrix */

  deter=fabs(matrix_det(a));
/*
  if (fabs(deter-1.0) > .0001)
    printf("\ndcm_to_rot:Warning: Invalid dir cos matrix (det=%f)\n",deter);
*/
/* determine direction of rotation */

  if (d2 == d1+1 || d2 == d1-2)
    rotdir = RT_FWD;
  else
    rotdir = RT_REV;

  d[0] = --d1;
  d[1] = --d2;
  d[2] = --d3;

/* determine sign change for nondiagonal elements */

  sign = type*rotdir;

/* 
  create new matrix with elements of a, rearranged to correspond with 
  desired type and rotation order (note - not a valid dircos matrix)
*/

  for (i=0; i<3; i++)
    {
      for (j=0; j<3; j++)
	{
	  if (type == RT_BODY)
	    {
	      bi = d[j];
	      bj = d[i];
	    }
	  else
	    {
	      bi = d[i];
	      bj = d[j];
	    }
	  b[i][j] = a[bi][bj];
	  if (i != j)
	    b[i][j] = sign*b[i][j];
	}
    }

/*
  use Kane's algorithm (Spacecraft Dynamics, p. 31) to determine
  rotation angles in the proper quadrant
*/

  if(fabs(fabs(b[2][0])-1.0) > 0.0001)
    {
      *r2 = asin(-b[2][0]);
      alpha = asin( b[2][1]/cos(*r2) );
      if (b[2][2] < 0.0)
	*r1 = PI - alpha;
      else
	*r1 = alpha;
      beta = asin( b[1][0]/cos(*r2) );
      if (b[0][0] < 0.0)
	*r3 = PI - beta;
      else
	*r3 = beta;
    }
  else
    {
      if (b[2][0] > 0.0)
	*r2 = -PI/2.;
      else
	*r2 = PI/2.;
      alpha = asin(-b[1][2]);
      if (b[1][1] < 0.0)
	*r1 = PI - alpha;
      else
	*r1 = alpha;
      *r3 = 0.0;
    }

}    



/***********************************************************************
NEW  dcm_to_rot_iji - takes a direction cos matrix a[3][3] and calculates
  body- or space-fixed rotation angles r1,r2,r3 in the order defined 
  by d1,d2,d3.
  Brings second body into orientation relative to the first by performing
  3 successive rotations which involve only *2* distinct unit vectors.
  For example, r1=x r2=y r3=x (only).

  SET UP NOW TO *ONLY* PERFORM RT_BODY-FIXED ROTATIONS! (wrote algorithm according
  to the space-two:1-2-1 (top p33 and matrix in appendix...) - Rearranges the 
  2-1-2 components to that of 1-2-1 for space-fixed and follows the given 
  algorithm -- USES RT_SPACE INSTEAD OF RT_BODY BECAUSE OF THE ROTATION ORDER - I 
  read that 3 body-fixed rotations should be multiplied in the reverse order, 
  but when I did this I got the matrix given in the book as space-fixed...

+++ so if you want body-fixed, use 1, not -1, and type in 1,2,1 to get y,x,y +++

  NOTE: this algorithm always assumes that  0 < r2 < PI, and
  calculates a solution accordingly.
*/

/*
#define RT_FWD  1
#define RT_REV -1
#define RT_SPACE 1
#define RT_BODY -1
*/

void dcm_to_rot_iji(double  a[3][3], 
                    int     type, 
                    int     d1, 
                    int     d2, 
                    int     d3, 
                    double* r1, 
                    double* r2, 
                    double* r3)

// double a[3][3];			/* input dir cos matrix */
// int type;			    /* 1 for space-fixed, -1 for body-fixed */
// int d1, d2, d3;			/* order for rotations */
// double *r1, *r2, *r3;    /* pointers to returned rotations */

{
  double b[3][3];		/* temp storage for rearranged el's of a */
  int d[3];			/* temp storage for --d1,--d2,--d3 */
  double sign;
  double alpha;
  double beta;
  double deter;
  double PI;
  // int i, j;
  // int bi, bj;
  int rotdir;

  PI = 4.0*atan(1.0);

/* check validity of input matrix */

  deter=fabs(matrix_det(a));
/*
  if (fabs(deter-1.0) > .0001)
    printf("\ndcm_to_rot:Warning: Invalid dir cos matrix (det=%f)\n",deter);
*/
/* determine direction of rotation */

  if (d2 == d1+1 || d2 == d1-2)
    rotdir = RT_FWD;
  else
    rotdir = RT_REV;

  d[0] = --d1;
  d[1] = --d2;
  d[2] = --d3;

/* determine sign change for nondiagonal elements */

  sign = type*rotdir;

/* 
  create new matrix with elements of a, rearranged to correspond with 
  desired type and rotation order (note - not a valid dircos matrix)
*/

  /*  for (i=0; i<3; i++)    */
  /*    {                    */
  /*      for (j=0; j<3; j++)*/
  /*	{                    */
  /*	  if (type == RT_BODY)  */
  /*	    {                */
  /*	      bi = d[j];     */
  /*	      bj = d[i];     */
  /*	    }                */
  /*	  else               */
  /*	    {                */
  /*	      bi = d[i];     */
  /*	      bj = d[j];     */
  /*	    }                */
  /*	  b[i][j] = a[bi][bj];      */
  /*	  if (i != j)               */
  /*	    b[i][j] = sign*b[i][j]; */
  /*	}                           */
  /*    }                           */


  /*move the components to the correct position for Kane's 1-2-1 algorithm*/
  b[0][0]=a[1][1];
  b[0][1]=a[1][0];
  b[0][2]=-a[1][2];
  b[1][0]=a[0][1];
  b[1][1]=a[0][0];
  b[1][2]=-a[0][2];
  b[2][0]=-a[2][1];
  b[2][1]=-a[2][0];
  b[2][2]=a[2][2];

/*
  use Kane's algorithm (Spacecraft Dynamics, p. 33) to determine
  rotation angles in the proper quadrant (for space-fixed 1-2-1)
*/

  if(fabs(fabs(b[0][0])-1.0) > 0.0001)
    {
      *r2 = acos(b[0][0]);
      alpha = asin(b[0][1]/sin(*r2) );
      if (b[0][2] < 0.0)
	*r1 = PI - alpha;
      else
	*r1 = alpha;
      beta = asin(b[1][0]/sin(*r2) );
      if (b[2][0] < 0.0)
	*r3 = beta;
      else
	*r3 = PI - beta;
    }
  else
    {
      if (b[0][0] > 0.0)
	*r2 = 0.0;
      else
	*r2 = PI;
      alpha = asin(-b[1][2]);
      if (b[1][1] < 0.0)
	*r1 = PI - alpha;
      else
	*r1 = alpha;
      *r3 = 0.0;
    }

}    


/***********************************************************************
  dcm_to_rot_jkj - takes a direction cos matrix a[3][3] and calculates
  body-fixed rotation angles r1,r2,r3.  Brings second body into orientation 
  relative to the first by performing 3 successive rotations which involve 
  only *2* distinct unit vectors.  For example, r1=Y r2=Z r3=Y (only).

  SET UP NOW TO *ONLY* PERFORM RT_BODY-FIXED ROTATIONS! (wrote algorithm according
  to the space-two:1-2-1 (top p33 and matrix in appendix...) - Rearranges the 
  2-3-2 components to that of 1-2-1 for space-fixed and follows the given 
  algorithm -- USES RT_SPACE INSTEAD OF RT_BODY BECAUSE OF THE ROTATION ORDER - I 
  read that 3 body-fixed rotations should be multiplied in the reverse order, 
  but when I did this I got the matrix given in the book as space-fixed...

  NOTE: this algorithm always assumes that  0 < r2 < PI, and
  calculates a solution accordingly.
*/

void dcm_to_rot_jkj(double  a[3][3], 
                    int     type, 
                    double* r1, 
                    double* r2, 
                    double* r3)

// double a[3][3];			/* input dir cos matrix */
// double *r1, *r2, *r3;		/* pointers to returned rotations */
// int type;                       /* 1 for space-fixed, -1 for body-fixed */

{
  double b[3][3];		/* temp storage for rearranged el's of a */
  // double sign;
  double alpha;
  double beta;
  double PI;
  // int i, j;
  // int bi, bj;
  // int rotdir;

  PI = 4.0*atan(1.0);

/* create new matrix with elements of a, rearranged to correspond with 
   desired type and rotation order (note - not a valid dircos matrix)*/
   /* move the components to the correct position for Kane's 1-2-1 algorithm*/
  
  b[0][0]=a[1][1];
  b[0][1]=a[1][2];
  b[0][2]=a[1][0];
  b[1][0]=a[2][1];
  b[1][1]=a[2][2];
  b[1][2]=a[2][0];
  b[2][0]=a[0][1];
  b[2][1]=a[0][2];
  b[2][2]=a[0][0];

/*
  use Kane's algorithm (Spacecraft Dynamics, p. 33) to determine
  rotation angles in the proper quadrant (for space-fixed 1-2-1)
*/

  if(type==1){                                 /* space-fixed */
    if(fabs(fabs(b[0][0])-1.0) > 0.0001){
      *r2 = acos(b[0][0]);
      alpha = asin(b[0][1]/sin(*r2) );
      if (b[0][2] < 0.0)
	*r1 = PI - alpha;
      else
	*r1 = alpha;
      beta = asin(b[1][0]/sin(*r2) );
      if (b[2][0] < 0.0)
	*r3 = beta;
      else
	*r3 = PI - beta;
    }
    else{
      if (b[0][0] > 0.0)
	*r2 = 0.0;
      else
	*r2 = PI;
      alpha = asin(-b[1][2]);
      if (b[1][1] < 0.0)
	*r1 = PI - alpha;
      else
	*r1 = alpha;
      *r3 = 0.0;
    }
  }
  else{                                        /* body-fixed */
    if(fabs(fabs(b[0][0])-1.0) > 0.0001){
      *r2 = acos(b[0][0]);
      alpha = asin(b[1][0]/sin(*r2) );
      if (b[2][0] < 0.0)
	*r1 = alpha;
      else
	*r1 = PI - alpha;
      beta = asin(b[0][1]/sin(*r2) );
      if (b[0][2] < 0.0)
	*r3 = PI - beta;
      else
	*r3 = beta;
    }
    else{
      if (b[0][0] > 0.0)
	*r2 = 0.0;
      else
	*r2 = PI;
      alpha = asin(b[2][1]);
      if (b[1][1] < 0.0)
	*r1 = PI - alpha;
      else
	*r1 = alpha;
      *r3 = 0.0;
    }
  }


}

/***********************************************************************
  rot_to_dcm - takes body- or space-fixed rotation angles r1,r2,r3
  (in the order defined by d1,d2,d3) and calculates a direction cos 
  matrix a[3][3]
*/

/*
#define RT_FWD  1
#define RT_REV -1
#define RT_SPACE 1
#define RT_BODY -1
*/

void rot_to_dcm(double  a[3][3], 
                int     type, 
                int     d1, 
                int     d2, 
                int     d3, 
                double  r1, 
                double  r2, 
                double  r3)

// double a[3][3];			/* output dir cos matrix */
// int type;			    /* 1 for space-fixed, -1 for body-fixed */
// int d1, d2, d3;			/* order for rotations */
// double r1, r2, r3;		/* input rotations */

{
  double b[3][3];		/* temp storage for rearranged el's of a */
  int d[3];			/* temp storage for --d1,--d2,--d3 */
  double sign;
  double deter;
  double c1, c2, c3;
  double s1, s2, s3;
  int i, j;
  int bi, bj;
  int rotdir;

/* determine direction of rotation */

  if (d2 == d1+1 || d2 == d1-2)
    rotdir = RT_FWD;
  else
    rotdir = RT_REV;

  d[0] = --d1;
  d[1] = --d2;
  d[2] = --d3;

/* determine sign change for nondiagonal elements */

  sign = type*rotdir;

/* calculate elements of dir cos matrix (assuming space3:123 rot) */

  c1 = cos(r1);
  c2 = cos(r2);
  c3 = cos(r3);
  s1 = sin(r1);
  s2 = sin(r2);
  s3 = sin(r3);

  b[0][0] =  c2*c3;
  b[0][1] =  s1*s2*c3 - sign*s3*c1;
  b[0][2] =  sign*c1*s2*c3 + s3*s1;
  b[1][0] =  sign*c2*s3;
  b[1][1] =  sign*s1*s2*s3 + c3*c1;
  b[1][2] =  c1*s2*s3 - sign*c3*s1;
  b[2][0] = -sign*s2;
  b[2][1] =  sign*s1*c2;
  b[2][2] =  c1*c2;

/* 
  create new matrix with elements of b, rearranged to correspond with 
  desired type and rotation order 
*/

  for (i=0; i<3; i++)
    {
      for (j=0; j<3; j++)
	{
	  if (type == RT_BODY)
	    {
	      bi = d[j];
	      bj = d[i];
	    }
	  else
	    {
	      bi = d[i];
	      bj = d[j];
	    }
	  a[bi][bj] = b[i][j];
	}
    }

/* check validity of ouput matrix */

  deter=fabs(matrix_det(a));
  if (fabs(deter-1.0) > .0001)
    printf("\nrot_to_dcm:Warning: Invalid dir cos matrix (det=%f)\n",deter);

}


/***************************************************************************

  orthonorm - uses Gram-Schmidt process to generate an orthonormal basis
  in R3 from 2 supplied non-colinear vectors.  Input vectors are supplied
  as pointers to VectorMT structures (defined at the top of this file).

  The vectors v1, v2 are processed into orthogonal unit vectors 
  u1,u2,u3 (which, together, make up the basis matrix) as follows:

  u1 is chosen as a unit vector in the direction of v1.
  u2 is chosen as a unit vector resulting from projecting v2 onto a plane
     prependicular to v1.
  u3 = u1 x u2.

  u1 is assigned to the axis designated by ax1 (1=x, 2=y, 3=z).
  u2 is assigned to the axis designated by ax2.

  

*/

void orthonorm(VectorMT*    v1, 
               VectorMT*    v2, 
               int          ax1, 
               int          ax2, 
               double       basis[3][3])

// VectorMT *v1, *v2;		/* input vector structures */
// int ax1, ax2;			/* axes corresponding to v1,v2 (1,2 or 3) */
// double basis[3][3];		/* output orthonormal basis matrix */

{

  int nrow, ncol;
  int axes[3];
  double temp;
  double tempmat[3][3];
  VectorMT u1, u2, u3;

/* u1 is unit vector in dir. of v1 */

  norm(v1,&u1);

/* use Gram-Schmidt to find u2 (from Intro. Linear Algebra, Kolman, p.161) */

  temp = dot(v2,&u1);
  u2.x = v2->x - temp*u1.x;
  u2.y = v2->y - temp*u1.y;
  u2.z = v2->z - temp*u1.z;
  norm(&u2,&u2);

/* determine order for cross-product resulting in R-H coord. system */

  if (ax2-ax1 == 1 || ax2-ax1 == -2)
    cross(&u1,&u2,&u3);
  else
    cross(&u2,&u1,&u3);

/* build basis (dir cos) matrix in correct order */

  build_matrix(&u1,&u2,&u3,tempmat);

  axes[0] = ax1-1;
  axes[1] = ax2-1;
  axes[2] = 5 - (ax1 + ax2);

  for (nrow=0; nrow<3; nrow++)
    for (ncol=0; ncol<3; ncol++)
      basis[nrow][axes[ncol]] = tempmat[nrow][ncol];


/* check validity of ouput matrix */

  temp=fabs(matrix_det(basis));
  if (fabs(temp-1.0) > .0001)
    printf("\northonorm:Warning: Invalid dir cos matrix (det=%f)\n",temp);

}

void build_cs(double    p1[3], 
              double    p2[3], 
              double    p3[3], 
              double    R[3][3])

// double p1[3], p2[3], p3[3], R[3][3];
/* make an orthogonal right handed CS with 
   unit vectors directed along each axis.  
   p1 is the origin, p2 is a point on the x-axis 
   and p3 is a point in the xy-plane.  
   R[xyz][0] is a unit vector directed along the x-axis
   R[xyz][1] is a unit vector directed along the y-axis
   R[xyz][2] is a unit vector directed along the z-axis
   So the unit vectors in the "R" matrix are in columns.
*/
{

  int xyz;
  double v1[3], v2[3], v3[3];
  VectorMT r1, r2, r3, r4;

  
  for(xyz=0;xyz<3;xyz++){
    v1[xyz] = p2[xyz] - p1[xyz]; 
    v2[xyz] = p3[xyz] - p1[xyz];
  }

  build_vector(v1,&r1);
  build_vector(v2,&r2);
  cross(&r1,&r2,&r3);
  cross(&r3,&r1,&r4);
  norm(&r1,&r1);
  norm(&r4,&r2);
  norm(&r3,&r3);
  strip_vector(&r1,v1);
  strip_vector(&r2,v2);
  strip_vector(&r3,v3);
  for(xyz=0;xyz<3;xyz++){
    R[xyz][0]=v1[xyz]; 
    R[xyz][1]=v2[xyz]; 
    R[xyz][2]=v3[xyz]; 
  }
}


void build_humerus_cs_r(double    p1[3], 
                        double    p2[3], 
                        double    p3[3], 
                        double    R[3][3])

// double p1[3], p2[3], p3[3], R[3][3];
/* make an orthogonal right handed CS with 
   unit vectors directed along each axis.  
   p1 is the origin, p2 is a point on the y-axis 
   and p3 is a point in the xy-plane.  For the RIGHT side.
   R[xyz][0] is a unit vector directed along the x-axis
   R[xyz][1] is a unit vector directed along the y-axis
   R[xyz][2] is a unit vector directed along the z-axis
   So the unit vectors in the "R" matrix are in columns.
*/
{

  int xyz;
  double v1[3], v2[3], v3[3];
  VectorMT r1, r2, r3, r4;

  
  for(xyz=0;xyz<3;xyz++){
    v1[xyz] = p2[xyz] - p1[xyz]; 
    v2[xyz] = p3[xyz] - p1[xyz];
  }

  build_vector(v1,&r1);
  build_vector(v2,&r2);
  cross(&r2,&r1,&r3);
  cross(&r3,&r1,&r4);
  norm(&r1,&r1);
  norm(&r4,&r2);
  norm(&r3,&r3);
  strip_vector(&r1,v1);
  strip_vector(&r2,v2);
  strip_vector(&r3,v3);
  for(xyz=0;xyz<3;xyz++){
    R[xyz][0]=v3[xyz]; 
    R[xyz][1]=v1[xyz]; 
    R[xyz][2]=v2[xyz]; 
  }
}



void build_humerus_cs_l(double    p1[3], 
                        double    p2[3], 
                        double    p3[3], 
                        double    R[3][3])

// double p1[3], p2[3], p3[3], R[3][3];
/* make an orthogonal right handed CS with 
   unit vectors directed along each axis.  
   p1 is the origin, p2 is a point on the y-axis 
   and p3 is a point in the xy-plane.  For the LEFT side
   R[xyz][0] is a unit vector directed along the x-axis
   R[xyz][1] is a unit vector directed along the y-axis
   R[xyz][2] is a unit vector directed along the z-axis
   So the unit vectors in the "R" matrix are in columns.
*/
{

  int xyz;
  double v1[3], v2[3], v3[3];
  VectorMT r1, r2, r3, r4;

  
  for(xyz=0;xyz<3;xyz++){
    v1[xyz] = p2[xyz] - p1[xyz]; 
    v2[xyz] = p3[xyz] - p1[xyz];
  }

  build_vector(v1,&r1);
  build_vector(v2,&r2);
  cross(&r1,&r2,&r3);
  cross(&r1,&r3,&r4);
  norm(&r1,&r1);
  norm(&r4,&r2);
  norm(&r3,&r3);
  strip_vector(&r1,v1);
  strip_vector(&r2,v2);
  strip_vector(&r3,v3);
  for(xyz=0;xyz<3;xyz++){
    R[xyz][0]=v2[xyz]; 
    R[xyz][1]=v1[xyz]; 
    R[xyz][2]=v3[xyz]; 
  }
}

/* */

void build_trunk_cs_r(  double    p1[3], 
                        double    p2[3], 
                        double    p3[3], 
                        double    R[3][3])

// double p1[3], p2[3], p3[3], R[3][3];
/* edited 092105 for new_eva_kinematics.c - for reflective 
   markers from markerless tracking subjects.
   make an orthogonal right handed CS with 
   unit vectors directed along each axis.  
   p2 is the origin, p1 is a point on the x-axis 
   and p3 is a point in the xz-plane.  For the RIGHT side. 
   R[xyz][0] is a unit vector directed along the x-axis
   R[xyz][1] is a unit vector directed along the y-axis
   R[xyz][2] is a unit vector directed along the z-axis
   So the unit vectors in the "R" matrix are in columns.
*/
{

  int xyz;
  double v1[3], v2[3], v3[3];
  VectorMT r1, r2, r3, r4;

  
  for(xyz=0;xyz<3;xyz++){
    v1[xyz] = p2[xyz] - p1[xyz]; 
    v2[xyz] = p2[xyz] - p3[xyz];
  }

  build_vector(v1,&r1);
  build_vector(v2,&r2);
  cross(&r1,&r2,&r3);
  cross(&r1,&r3,&r4);
  norm(&r1,&r1);
  norm(&r4,&r2);
  norm(&r3,&r3);
  strip_vector(&r1,v1);
  strip_vector(&r2,v2);
  strip_vector(&r3,v3);
  for(xyz=0;xyz<3;xyz++){
    R[xyz][0]=v1[xyz]; 
    R[xyz][1]=v3[xyz]; 
    R[xyz][2]=v2[xyz]; 
  }
}

void build_trunk_cs_l(  double    p1[3], 
                        double    p2[3], 
                        double    p3[3], 
                        double    R[3][3])

// double p1[3], p2[3], p3[3], R[3][3];
/* edited 092105 for new_eva_kinematics.c - for reflective 
   markers from markerless tracking subjects.
   make an orthogonal right handed CS with 
   unit vectors directed along each axis.  
   p2 is the origin, p1 is a point on the x-axis 
   and p3 is a point in the xz-plane.  For the LEFT side,
   but gives a flip around the y-axis from the right side
   coord system created above.  Need to make z values negative
   in output of main program so that z means "out the
   back" for both left and right.
   R[xyz][0] is a unit vector directed along the x-axis
   R[xyz][1] is a unit vector directed along the y-axis
   R[xyz][2] is a unit vector directed along the z-axis
   So the unit vectors in the "R" matrix are in columns.
*/
{

  int xyz;
  double v1[3], v2[3], v3[3];
  VectorMT r1, r2, r3, r4;

  
  for(xyz=0;xyz<3;xyz++){
    v1[xyz] = p2[xyz] - p1[xyz]; 
    v2[xyz] = p2[xyz] - p3[xyz];
  }

  build_vector(v1,&r1);
  build_vector(v2,&r2);
  cross(&r2,&r1,&r3);
  cross(&r1,&r3,&r4);
  norm(&r1,&r1);
  norm(&r4,&r2);
  norm(&r3,&r3);
  strip_vector(&r1,v1);
  strip_vector(&r2,v2);
  strip_vector(&r3,v3);
  for(xyz=0;xyz<3;xyz++){
    R[xyz][0]=v1[xyz]; 
    R[xyz][1]=v3[xyz]; 
    R[xyz][2]=v2[xyz]; 
  }
}

void build_vector_cs(   double    p1[3], 
                        double    p2[3], 
                        double    p3[3], 
                        double    p4[3], 
                        double    R[3][3])

// double p1[3], p2[3], p3[3], p4[3];
// double R[3][3];
/* build right handed CS using 4 points as input*/

{
  int xyz;
  double v1[3], v2[3], v3[3];
  VectorMT r1, r2, r3, r4, r5, r6, r7, r8;

  build_vector(p1, &r1);
  build_vector(p2, &r2);
  build_vector(p3, &r3);
  build_vector(p4, &r4);
  sub_vector(&r3,&r1,&r5);
  sub_vector(&r4,&r2,&r6);
  cross(&r5,&r6,&r7);
  cross(&r7,&r5,&r8);
  norm(&r5,&r5);
  norm(&r8,&r8);
  norm(&r7,&r7);
  strip_vector(&r5,v1);
  strip_vector(&r8,v2);
  strip_vector(&r7,v3);
  for(xyz=0;xyz<3;xyz++){
    R[xyz][0]=v1[xyz]; 
    R[xyz][1]=v2[xyz]; 
    R[xyz][2]=v3[xyz]; 
  }
}


void build_vector_cs_hum(   double    p1[3], 
                            double    p2[3], 
                            double    p3[3], 
                            double    R[3][3])

// double p1[3], p2[3], p3[3];
// double R[3][3];
/* build right handed CS using 3 points as input.
   It has been edited 092105 to work in 
   new_eva_kinematics.c for markerless tracking 
   subjects' reflective marker data, on both sides,
   left and right. p1 is origin, p3 is on xy-plane */

{
  int xyz;
  double v1[3], v2[3], v3[3];
  VectorMT r1, r2, r3, /*r4,*/ r5, r6, r7, r8;

  build_vector(p1, &r1);
  build_vector(p2, &r2);
  build_vector(p3, &r3);
  sub_vector(&r3,&r1,&r5);
  sub_vector(&r2,&r1,&r6);
  cross(&r6,&r5,&r7);
  cross(&r5,&r7,&r8);
  norm(&r5,&r5);
  norm(&r8,&r8);
  norm(&r7,&r7);
  strip_vector(&r5,v1);
  strip_vector(&r8,v2);
  strip_vector(&r7,v3);
  for(xyz=0;xyz<3;xyz++){
    R[xyz][0]=v2[xyz]; 
    R[xyz][1]=v1[xyz]; 
    R[xyz][2]=v3[xyz]; 
  }
}


void build_vector_cs_hum_RC(double    p1[3], 
                            double    p2[3], 
                            double    p3[3], 
                            double    p4[3], 
                            double    R[3][3])

// double p1[3], p2[3], p3[3], p4[3];
// double R[3][3];
/* build right handed CS using 4 points as input.
   It has been edited 100505 to work in 
   kinematics_RCrepair_humnotreltofirstframe.c for 
   markerless tracking subjects' anatomic marker 
   data, on both sides. */

{
  int xyz;
  double v1[3], v2[3], v3[3];
  VectorMT r1, r2, r3, r4, r5, r6, r7, r8;

  build_vector(p1, &r1);
  build_vector(p2, &r2);
  build_vector(p3, &r3);
  build_vector(p4, &r4);
  sub_vector(&r3,&r1,&r5);
  sub_vector(&r4,&r2,&r6);
  cross(&r6,&r5,&r7);
  cross(&r5,&r7,&r8);
  norm(&r5,&r5);
  norm(&r8,&r8);
  norm(&r7,&r7);
  strip_vector(&r5,v1);
  strip_vector(&r8,v2);
  strip_vector(&r7,v3);
  for(xyz=0;xyz<3;xyz++){
    R[xyz][0]=v2[xyz]; 
    R[xyz][1]=v1[xyz]; 
    R[xyz][2]=v3[xyz]; 
  }
}


void build_vertebrae_cs(    double    p1[3], 
                            double    p2[3], 
                            double    p3[3], 
                            double    p4[3], 
                            double    R[3][3])

// double p1[3], p2[3], p3[3], p4[3];
// double R[3][3];
/* build right handed CS using 4 points as input */
/* p1 and p2 medial and lateral edges of the body, p3 and p4 are anterior and posterior */

{
  int xyz;
  double v1[3], v2[3], v3[3];                   //, v4[3];
  VectorMT r1, r2, r3, r4, /*r5, r6,*/ r7, r8;  //, r9;


  build_vector(p1, &r1);
  build_vector(p2, &r2);
  build_vector(p3, &r3);
  build_vector(p4, &r4);
  sub_vector(&r2,&r1,&r7);
  sub_vector(&r3,&r4,&r8);
  cross(&r8,&r7,&r1);
  cross(&r7,&r1,&r2);
  norm(&r1,&r1);
  norm(&r2,&r2);
  norm(&r7,&r7);
  strip_vector(&r2,v1);
  strip_vector(&r7,v2);
  strip_vector(&r1,v3);
  for(xyz=0;xyz<3;xyz++){
    R[xyz][0]=v1[xyz]; 
    R[xyz][1]=v2[xyz]; 
    R[xyz][2]=v3[xyz]; 
  }
}

void build_shoulder_cs( double    p1[3], 
                        double    p2[3], 
                        double    p3[3], 
                        double    R[3][3])

// double p1[3], p2[3], p3[3], R[3][3];
/* make an orthogonal right handed CS with 
   unit vectors directed along each axis.  
   Follow the method from van der Helm.
   R[xyz][0] is a unit vector directed along the x-axis
   R[xyz][1] is a unit vector directed along the y-axis
   R[xyz][2] is a unit vector directed along the z-axis
   So the unit vectors in the "R" matrix are in columns.
   p1 is on x-axis, p2 is origin, p3 is in xy plane.
*/
{

  int xyz;
  double v1[3], v2[3], v3[3];
  VectorMT r1, r2, r3, r4;

  
  for(xyz=0;xyz<3;xyz++){
    v1[xyz] = p2[xyz] - p1[xyz]; 
    v2[xyz] = p2[xyz] - p3[xyz];
  }

  build_vector(v1,&r1);
  build_vector(v2,&r2);
  cross(&r1,&r2,&r3);
  cross(&r3,&r1,&r4);
  norm(&r1,&r1);
  norm(&r4,&r2);
  norm(&r3,&r3);
  strip_vector(&r1,v1);
  strip_vector(&r2,v2);
  strip_vector(&r3,v3);
  for(xyz=0;xyz<3;xyz++){
    R[xyz][0]=v1[xyz]; 
    R[xyz][1]=v2[xyz]; 
    R[xyz][2]=v3[xyz]; 
  }
}

/* */

void projection_angles(double   R1[3][3], 
                       double   R2[3][3], 
                       double*  ang_x, 
                       double*  ang_y, 
                       double*  ang_z)
     /* subroutine to find projection angles between two coordinate systems */
     /* calculates orientation of "final" dir cos matrix (R2 in ihat, jhat, khat) */
     /* rel to "intial" dir cos matrix (R1 in X,Y,Z)*/

// double R1[3][3], R2[3][3];       /* R1 is "initial" dir cos matrix; R2 is "final" dir cos matrix */
// double *ang_x, *ang_y, *ang_z;   /* pointers to returned rotations */


{
  double R3[3][3], R[3][3];
  double PI= 3.1415926535;
  transpose(R1,R3);
  mult_matrix(R3,R2,R);
  *ang_x=(180./PI)*atan(R[2][1]/R[1][1]);
  *ang_y=-(180./PI)*atan(R[2][0]/R[0][0]);
  *ang_z=(180./PI)*atan(R[1][0]/R[0][0]);

}

void projection_angles2(double   R1[3][3], 
                        double   R2[3][3], 
                        double*  ang_x, 
                        double*  ang_y, 
                        double*  ang_z)

     /* subroutine to find projection angles between two coordinate systems */
     /* calculates orientation of "final" dir cos matrix (R2 in ihat, jhat, khat) */
     /* rel to "intial" dir cos matrix (R1 in X,Y,Z)*/

// double R1[3][3], R2[3][3];  /* R1 is "initial" dir cos matrix; R2 is "final" dir cos matrix */
// double *ang_x, *ang_y, *ang_z;		/* pointers to returned rotations */


{
  double v1, v2;    //, v3;

  double R3[3][3], R[3][3];
  double PI= 3.1415926535;

  transpose(R1,R3);
  mult_matrix(R3,R2,R);
  

  v1=sqrt(R[0][1]*R[0][1]+R[1][1]*R[1][1]);
  v2=sqrt(R[0][0]*R[0][0]+R[1][0]*R[1][0]);
    

  *ang_x=(180./PI)*atan(R[2][1]/v1);
  *ang_y=-(180./PI)*atan(R[2][0]/v2);
  *ang_z=-(180./PI)*asin(R[0][1]/v1);

  
}

void reorder(double p1[3], 
             double p2[3], 
             double p3[3], 
             int*   i, 
             int*   j, 
             int*   k)

     /* subroutine ro reorder vertices so that longest vector is from 1 to 2 */
     /* and second longest vector is from 1 to 3 */

// double p1[3], p2[3], p3[3];
// int *i, *j, *k;

{
  double d1, d2, d3;

  d1=distance((VectorMT*)p1, (VectorMT*)p2);    // NOTE: (VectorMT*)p1 is synonymous with double p1[3]
  d2=distance((VectorMT*)p1, (VectorMT*)p3);
  d3=distance((VectorMT*)p2, (VectorMT*)p3);

  if(d1>d2&&d2>d3){
    *i=0;
    *j=1;
    *k=2;
  }
  if(d1>d3&&d3>d2){
    *i=1;
    *j=0;
    *k=2;
  }
  if(d3>d2&&d2>d1){
    *i=2;
    *j=1;
    *k=0;
  }
  if(d3>d1&&d1>d2){
    *i=2;
    *j=0;
    *k=1;
  }
  if(d2>d1&&d1>d3){
    *i=0;
    *j=2;
    *k=1;
  }
  if(d2>d3&&d3>d1){
    *i=1;
    *j=2;
    *k=0;
  }


}
