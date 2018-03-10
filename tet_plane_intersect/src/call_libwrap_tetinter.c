#include <stdlib.h>
#include <stdio.h>
#include "plane_normal_tetrahedron_intersect.h"


int main(int argc, char **argv) {

  // Define plane and tetrahedra
  double pp[3] = {0.5, 0.0, 0.0};
  double normal[3] = {1.0, 0.0, 0.0};
  int n_tets = 3;
  double tets[36] = \
      { 0.0696218 ,  0.51675287,  0.27622866,
        0.77795185,  0.35896   ,  0.76889514,
        0.93955095,  0.46952051,  0.50414596,
        0.88091001,  0.0994195 ,  0.58077994,
        0.22094245,  0.2076138 ,  0.97860943,
        0.07205741,  0.20692974,  0.03908092,
        0.42881122,  0.14083737,  0.73928276,
        0.81604653,  0.96100499,  0.73425617,
        0.36762177,  0.19788064,  0.32692846,
        0.74746389,  0.26432934,  0.82284067,
        0.59302365,  0.03474258,  0.48759837,
        0.90889265,  0.50532202,  0.09855027 };

  // Prealloc output arrays
  double *pint = (double *) malloc(n_tets * 18 * sizeof(double));
  double *int_area = (double *) malloc(n_tets * 6 * sizeof(double));
  double *int_center = (double *) malloc(n_tets * 6 * sizeof(double));
  int *tri_num = (int *) malloc(n_tets * sizeof(int));
  int *tri_ptr = (int *) malloc((n_tets + 1) * sizeof(int));

  // Call library function for batch intersections
  batch_intersect_tris(pp, normal, tets, n_tets, pint, int_area, int_center,\
      tri_num, tri_ptr);

  int n_tris = tri_ptr[n_tets];
  for (int i = 0; i < n_tris; i++) {
    printf("t %d, area %1.3e, center {%+1.3e %+1.3e %+1.3e}\n", i, int_area[i],\
        int_center[3*i], int_center[3*i+1], int_center[3*i+2]);
    for (int j = 0; j < 3; j++)
      printf("  -> p %d {%+1.3e %+1.3e %+1.3e}\n", j, pint[9*i+3*j],\
          pint[9*i+3*j+1], pint[9*i+3*j+2]);
  }

  // Free heap arrays
  free(pint);
  free(int_area);
  free(int_center);
  free(tri_num);
  free(tri_ptr);

  return 0;

}
