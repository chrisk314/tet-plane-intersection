#ifndef PLANE_NORMAL_TETRAHEDRON_INTERSECT_H_
#define PLANE_NORMAL_TETRAHEDRON_INTERSECT_H_

#define TOL 1.0e-15

double parallelogram_area_3d(double (*)[3]);
double quad_area_3d(double (*)[3]);
double tri_area_3d(double (*)[3]);
void set_inter_poly_center(double (*)[3], int, double *);
void plane_normal_tetrahedron_intersect (double [3], double [3], double (*)[3],\
    double (*)[3], int *);
void batch_intersect_polys(double [3], double [3], double *, int , double *,\
    double *, double *, int *, int *);
void batch_intersect_tris(double [3], double [3], double *, int, double *,\
    double *, double *, int *, int *);
void batch_intersect_tris_omp(double [3], double [3], double *, int, double *,\
    double *, double *, int *, int *);

#endif // PLANE_NORMAL_TETRAHEDRON_INTERSECT_H_
