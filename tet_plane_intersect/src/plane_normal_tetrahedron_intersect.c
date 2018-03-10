/*
 * This code is adapted from plane_normal_tetrahedron_intersect.m by John Burkardt.
 * Ported to C with Python wrapping by Chris Knight on 09 March 2018.
 *
 * The functions plane_normal_tetrahedron_intersect, parallelogram_area_3d, and
 * quad_area_3d are direct ports of John Burkardt's Matlab code. All the rest is my
 * original work.
 *
 * The derived functions batch_intersect_tris, and batch_intersect_polys can be used 
 * to find the intersections, as a set of triangles or polygons respectively, for a 
 * list of tetrahedra in batch.
 *
 * The original code can be found at this web page at the time of writing:
 * https://people.sc.fsu.edu/~jburkardt/m_src/tetrahedron_slice_animate/tetrahedron_slice_animate.html
 *
 * Original description including licence from John Burkardt's code is
 * reproduced below (the below description contains Matlab syntax/notation):
 * 
 * %%  PLANE_NORMAL_TETRAHEDRON_INTERSECT intersects a plane and a tetrahedron.
 * %%
 * %%  Discussion:
 * %%
 * %%    The intersection of a plane and a tetrahedron is one of:
 * %%    0) empty
 * %%    1) a single point
 * %%    2) a single line segment
 * %%    3) a triangle
 * %%    4) a quadrilateral.
 * %%
 * %%    In each case, the region of intersection can be described by the
 * %%    corresponding number of points.  In particular, cases 2, 3 and 4
 * %%    are described by the vertices that bound the line segment, triangle,
 * %%    or quadrilateral.
 * %%
 * %%    The normal form of a plane is:
 * %%
 * %%      PP is a point on the plane,
 * %%      N is a normal vector to the plane.
 * %%
 * %%    The form of a tetrahedron is
 * %%
 * %%      T(1:3,1:4) contains the coordinates of the vertices.
 * %%
 * %%  Licensing:
 * %%
 * %%    This code is distributed under the GNU LGPL license.
 * %%
 * %%  Modified:
 * %%
 * %%    24 June 2010
 * %%
 * %%  Author:
 * %%
 * %%    John Burkardt
 * %%
 * %%  Parameters:
 * %%
 * %%    Input, real PP(3), a point on the plane.
 * %%
 * %%    Input, real NORMAL(3), a normal vector to the plane.
 * %%
 * %%    Input, real T(3,4), the tetrahedron vertices.
 * %%
 * %%    Output, integer INT_NUM, the number of intersection
 * %%    points returned.  This will be 0, 1, 2, 3 or 4.
 * %%
 * %%    Output, real PINT(3,4), the coordinates of the
 * %%    intersection points.
 */

#include <string.h>
#include <math.h>
#include "plane_normal_tetrahedron_intersect.h"


double parallelogram_area_3d(double (*p)[3])
/*
 *  Parameters:
 *
 *    double q[4][3], parallelogram vertex coordinates
 *
 *  Returns:
 *    
 *    double parallelogram area
 */
{
  double cross[3];
  cross[0] = (p[1][1] - p[0][1]) * (p[2][2] - p[0][2]) \
           - (p[1][2] - p[0][2]) * (p[2][1] - p[0][1]);
                                                     
  cross[1] = (p[1][2] - p[0][2]) * (p[2][0] - p[0][0]) \
           - (p[1][0] - p[0][0]) * (p[2][2] - p[0][2]);
                                                     
  cross[2] = (p[1][0] - p[0][0]) * (p[2][1] - p[0][1]) \
           - (p[1][1] - p[0][1]) * (p[2][0] - p[0][0]);

  double mag = 0.;
  for (int i = 0; i < 3; i++)
    mag += cross[i] * cross[i];

  return sqrt(mag);
}


double quad_area_3d(double (*q)[3]) 
/*
 *  Parameters:
 *
 *    double q[4][3], quadrilateral vertex coordinates
 *
 *  Returns:
 *    
 *    double quadrilateral area
 */
{
  double p[4][3];
  for (int i = 0; i < 3; i++){
    p[0][i] = 0.5 * (q[0][i] + q[1][i]);
    p[1][i] = 0.5 * (q[1][i] + q[2][i]);
    p[2][i] = 0.5 * (q[2][i] + q[3][i]);
    p[3][i] = 0.5 * (q[3][i] + q[0][i]);
  }
  return 2.0 * parallelogram_area_3d(p);
}


double tri_area_3d(double (*t)[3]) 
/*
 *  Parameters:
 *
 *    double t[3][3], triangle vertex coordinates
 *
 *  Returns:
 *    
 *    double triangle area
 */
{
  return 0.5 * parallelogram_area_3d(t);
}


void set_inter_poly_center(double (*p)[3], int int_num, double *center)
/*
 *  Parameters:
 *
 *    double p[int_num][3], intersection vertex coordinates
 *
 *    int int_num, number of intersection points
 *
 *    double center[3], center of mass of intersection points
 */
{
  double inv_int_num = 1. / int_num;
  for (int i = 0; i < 3; i++) {
    center[i] = 0.;
    for (int j = 0; j < int_num; j++) {
      center[i] += p[j][i];
    }
    center[i] *= inv_int_num;
  }
  return;
}


void plane_normal_tetrahedron_intersect(double pp[3], double normal[3],\
    double (*t)[3], double (*pint)[3], int *int_num)
/*
 *  Parameters:
 *
 *    double pp[3], a point on the plane.
 *
 *    double normal[3], a normal vector to the plane.
 *
 *    double t[4][3], the tetrahedron vertices.
 *
 *    integer int_num, the number of intersection points returned.  This will be
 *    0, 1, 2, 3 or 4.
 *
 *    double pint[4][3], the coordinates of the intersection points.
 */
{
  *int_num = 0;
  memset(pint, 0., 12 * sizeof(double));

  // Normalise input normal vector
  double mag = 0.;
  for (int i = 0; i < 3; i++)
    mag += normal[i] * normal[i];
  mag = sqrt(mag);
  for (int i = 0; i < 3; i++)
    normal[i] /= mag;

  // Calculate distance from origin to the plane.
  double dpp = 0.;
  for (int i = 0; i < 3; i++)
    dpp += normal[i] * pp[i];

  // d[i] is positive, zero, or negative if vertex i is above, on, or below the
  // plane.
  double d[4];
  for (int i = 0; i < 4; i++){
    d[i] = 0.;
    for (int j = 0; j < 3; j++)
      d[i] += normal[j] * t[i][j];
    d[i] -= dpp;
  }

  // If all d are positive or negative, no intersection.
  if (
      (d[0] < 0. && d[1] < 0. && d[2] < 0 && d[3] < 0.) || 
      (d[0] > 0. && d[1] > 0. && d[2] > 0 && d[3] > 0.)
     )
  {
    return;
  }

  // Points with zero distance are automatically added to the list.
  //
  // For each point with nonzero distance, seek another point with opposite sign
  // and higher index, and compute the intersection of the line between those
  // points and the plane.
  for (int i = 0; i < 4; i++) {
    if ( fabs(d[i]) < TOL ) {
      memcpy(pint[(*int_num)], t[i], 3 * sizeof(double));
      (*int_num)++;
    }
    else{
      for (int j = i + 1; j < 4; j++) {
        if (d[i] * d[j] < 0.) {
          double inv_dist = 1. / (d[i] - d[j]);
          for (int k = 0; k < 3; k++) {
            pint[(*int_num)][k] = (d[i] * t[j][k] - d[j] * t[i][k]) * inv_dist;;
          }
          (*int_num)++;
        }
      }
    }
  }

  //  If four points were found, order them properly.
  if ( (*int_num) == 4 ) {
    double area1 = quad_area_3d(pint);
    double pint2[4][3];
    memcpy(pint2, pint, 6 * sizeof(double));
    memcpy(pint2[2], pint[3], 3 * sizeof(double));
    memcpy(pint2[3], pint[2], 3 * sizeof(double));
    double area2 = quad_area_3d(pint2);
    if (area1 < area2)
      memcpy(pint, pint2, 12 * sizeof(double));
  }

  return;
}


void batch_intersect_polys(double pp[3], double normal[3], double *T_list, int n_T,\
    double *pint, double *int_area, double *int_center, int *int_num, int *int_ptr)
/*
 *  Parameters:
 *
 *    double pp[3], a point on the plane.
 *
 *    double normal[3], a normal vector to the plane.
 *
 *    double *Tet_list[n_Tet][4][3], ptr to list of vertices for n_Tet tetrahedra.
 *
 *    int n_Tet, number of tetrahedra in list Tet_list.
 *
 *    double pint[n_Tet][4][3], the coordinates of the intersection points.
 *    There can be up to 4 intersection points per tetrahedron.
 *
 *    double int_area[n_Tet], the areas of the intersection polygons.
 *
 *    double int_center[n_Tet][3], the center coordinates of the intersection
 *    polygon.
 *
 *    integer int_num[n_Tet], the number of intersection points for each
 *    tetrahedron.  This will be 0, 1, 2, 3, or 4.
 *
 *    integer int_ptr[n_Tet+1], offsets of the starting point of the intersection
 *    point data for each tetrahedron.
 */
{
  int count = 0;

  int_ptr[0] = 0;

  for (int i = 0; i < n_T; i++){
    plane_normal_tetrahedron_intersect(pp, normal, (double (*)[3]) &T_list[12*i],\
        (double (*)[3]) &pint[3*int_ptr[i]], &int_num[i]);
    int_ptr[i+1] = int_ptr[i] + int_num[i];

    // Set intersection polygon areas and centers
    if (int_num[i] == 4) {
      int_area[count] = quad_area_3d((double (*)[3]) &pint[3*int_ptr[i]]);
      set_inter_poly_center((double (*)[3]) &pint[3*int_ptr[i]], int_num[i],\
          &int_center[3*count]);
      count++;
    }
    else if (int_num[i] == 3) {
      int_area[count] = tri_area_3d((double (*)[3]) &pint[3*int_ptr[i]]);
      set_inter_poly_center((double (*)[3]) &pint[3*int_ptr[i]], int_num[i],\
          &int_center[3*count]);
      count++;
    }
  }
  return;
}


void batch_intersect_tris(double pp[3], double normal[3], double *Tet_list, int n_Tet,\
    double *pint, double *int_area, double *int_center, int *tri_num, int *tri_ptr)
/*
 *  Parameters:
 *
 *    double pp[3], a point on the plane.
 *
 *    double normal[3], a normal vector to the plane.
 *
 *    double *Tet_list[n_Tet][4][3], ptr to list of vertices for n_Tet tetrahedra.
 *
 *    int n_Tet, number of tetrahedra in list Tet_list.
 *
 *    double pint[2*n_Tet][3][3], the coordinates of the intersection points.
 *    There can be up to 2 intersection triangles per tetrahedron.
 *
 *    double int_area[2*n_Tet], the areas of the intersection triangles.
 *
 *    double int_center[2*n_Tet][3], the center coordinates of the intersection
 *    triangles.
 *
 *    integer tri_num[n_Tet], the number of intersection triangle for each
 *    tetrahedron.  This will be 0, 1, or 2.
 *
 *    integer tri_ptr[n_Tet+1], offsets of the starting point of the intersection
 *    point data for each tetrahedron.
 */
{
  int int_num;
  int tri_count = 0;

  tri_ptr[0] = 0;
  
  for (int i = 0; i < n_Tet; i++){
    tri_num[i] = 0;

    plane_normal_tetrahedron_intersect(pp, normal, (double (*)[3]) &Tet_list[12*i],\
        (double (*)[3]) &pint[9*tri_count], &int_num);

    // Set intersection triangle areas and centers
    if (int_num == 4) {
      // Split quad into two triangles
      memcpy(&pint[9*(tri_count+1)+3], &pint[9*tri_count], 3 * sizeof(double));
      memcpy(&pint[9*(tri_count+1)+6], &pint[9*tri_count+6], 3 * sizeof(double));
      int_area[tri_count] = tri_area_3d((double (*)[3]) &pint[9*tri_count]);
      set_inter_poly_center((double (*)[3]) &pint[9*tri_count], 3,\
          &int_center[3*tri_count]);
      tri_count++;
      int_area[tri_count] = tri_area_3d((double (*)[3]) &pint[9*tri_count]);
      set_inter_poly_center((double (*)[3]) &pint[9*tri_count], 3,\
          &int_center[3*tri_count]);
      tri_count++;
      tri_num[i] = 2;
    }
    else if (int_num == 3) {
      int_area[tri_count] = tri_area_3d((double (*)[3]) &pint[9*tri_count]);
      set_inter_poly_center((double (*)[3]) &pint[9*tri_count], 3,\
          &int_center[3*tri_count]);
      tri_count++;
      tri_num[i] = 1;
    }
    tri_ptr[i+1] = tri_ptr[i] + tri_num[i];
  }
  return;
}
