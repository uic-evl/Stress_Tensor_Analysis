// defs for MathTools.cc & RotTools.cc

#pragma once

typedef struct VectorMT
{		         /* x,y,z vector coords */
  double x;
  double y;
  double z;
} VectorMT;

void norm(VectorMT* v, VectorMT* normv);
void cross(VectorMT* v1, VectorMT* v2, VectorMT* vc);
void copy_vector(VectorMT* v1, VectorMT* vout);
void add_vector(VectorMT* v1, VectorMT* v2, VectorMT* vout);
void sub_vector(VectorMT* v1, VectorMT* v2, VectorMT* vout);
void build_matrix(VectorMT *v1, VectorMT *v2, VectorMT *v3, double m[3][3]);
void build_matrix2(VectorMT *v1, VectorMT *v2, VectorMT *v3, double m[3][3]);
void strip_matrix(VectorMT *v1, VectorMT *v2, VectorMT *v3, double m[3][3]);
void strip_matrix2(VectorMT *v1, VectorMT *v2, VectorMT *v3, double m[3][3]);
void mult_matrix(double m1[3][3], double m2[3][3], double mprod[3][3]);
void build_vector(double v[3], VectorMT *vout);
void strip_vector(VectorMT *v, double vout[3]);
void matxvec(double m[3][3], VectorMT *v, VectorMT *vout);
void vecxmat(double m[3][3], VectorMT *v, VectorMT *vout);
void vecxmat2(double m[3][3], VectorMT *v, VectorMT *vout);
void transpose(double m[3][3], double mout[3][3]);
void copy_matrix(double m1[3][3], double mout[3][3]);
void vproj(VectorMT* v1, VectorMT* v2, VectorMT* v3);
double dot(VectorMT* v1, VectorMT* v2);
double distance(VectorMT* v1, VectorMT* v2);
double matrix_det(double a[3][3]);
double vmag(VectorMT* v);
double vangle(VectorMT* v1, VectorMT* v2);
void find_centroid(double v1[3], double v2[3], double v3[3], double v4[3]);
void find_triangle_area(double v1[3], double v2[3], double v3[3], double area);
void print_matrix(double m[3][3]);
double vproj_value(VectorMT* v1, VectorMT* v2);
void sxvector(double s, VectorMT* v, VectorMT* vout);
void convert_vector(double v1[3], double v2[3], double m[3][3], double vout[3]);

#define RT_FWD      1
#define RT_REV     -1
#define RT_SPACE    1
#define RT_BODY    -1

void dcm_to_rot(double  a[3][3], 
                int     type, 
                int     d1, 
                int     d2, 
                int     d3, 
                double* r1, 
                double* r2, 
                double* r3);

void dcm_to_rot_iji(double  a[3][3], 
                    int     type, 
                    int     d1, 
                    int     d2, 
                    int     d3, 
                    double* r1, 
                    double* r2, 
                    double* r3);

void dcm_to_rot_jkj(double  a[3][3], 
                    int     type, 
                    double* r1, 
                    double* r2, 
                    double* r3);

void rot_to_dcm(double  a[3][3], 
                int     type, 
                int     d1, 
                int     d2, 
                int     d3, 
                double  r1, 
                double  r2, 
                double  r3);

void orthonorm(VectorMT*    v1, 
               VectorMT*    v2, 
               int          ax1, 
               int          ax2, 
               double       basis[3][3]);

void build_cs(double    p1[3], 
              double    p2[3], 
              double    p3[3], 
              double    R[3][3]);

void build_humerus_cs_r(double    p1[3], 
                        double    p2[3], 
                        double    p3[3], 
                        double    R[3][3]);

void build_humerus_cs_l(double    p1[3], 
                        double    p2[3], 
                        double    p3[3], 
                        double    R[3][3]);

void build_trunk_cs_r(  double    p1[3], 
                        double    p2[3], 
                        double    p3[3], 
                        double    R[3][3]);

void build_trunk_cs_l(  double    p1[3], 
                        double    p2[3], 
                        double    p3[3], 
                        double    R[3][3]);

void build_vector_cs(   double    p1[3], 
                        double    p2[3], 
                        double    p3[3], 
                        double    p4[3], 
                        double    R[3][3]);

void build_vector_cs_hum(   double    p1[3], 
                            double    p2[3], 
                            double    p3[3], 
                            double    R[3][3]);

void build_vector_cs_hum_RC(double    p1[3], 
                            double    p2[3], 
                            double    p3[3], 
                            double    p4[3], 
                            double    R[3][3]);

void build_vertebrae_cs(    double    p1[3], 
                            double    p2[3], 
                            double    p3[3], 
                            double    p4[3], 
                            double    R[3][3]);

void build_shoulder_cs( double    p1[3], 
                        double    p2[3], 
                        double    p3[3], 
                        double    R[3][3]);

void projection_angles(double   R1[3][3], 
                       double   R2[3][3], 
                       double*  ang_x, 
                       double*  ang_y, 
                       double*  ang_z);

void projection_angles2(double   R1[3][3], 
                        double   R2[3][3], 
                        double*  ang_x, 
                        double*  ang_y, 
                        double*  ang_z);

void reorder(double p1[3], 
             double p2[3], 
             double p3[3], 
             int*   i, 
             int*   j, 
             int*   k);

