# Tet Plane Intersect
Python wrapper around a C library for calculating intersections between tetrahedra and planes. This
code is based on John Burkardt's [Matlab
code](https://people.sc.fsu.edu/~jburkardt/m_src/tetrahedron_slice_animate/tetrahedron_slice_animate.html) which calculates the
intersection points between a single tetrahedron and a plane. This functionality is extended to
enable calculation of intersections for multiple tetrahedra with a given plane in batch.

### Installation
To install simply run the below commands

```
git clone git@github.com:chrisk314/tet-plane-intersection.git
cd tet-plane-intersect
pip install .
```

### Calculating intersections
The function `tet_plane_inter_tris_batch` can be used to calculate the intersections between 
a collection of tetrahedra and a plane in batch. The intersection of a tetrahedron with a plane 
has 0, 1, 2, 3, or 4 points. This function converts any 4-point intersections into two 3-point
intersections, i.e., two triangles. This provides the convenience of working with the intersections
in the form of triangles only.

The below example code calls `tet_plane_inter_tris_batch` to obtain the triangles
produced by the intersection between the tetrahedra geometry specified in the `tetrahedra` numpy
array and the plane specified by the point `pp` and normal `normal`.

```python
import numpy as np
from tet_plane_intersect import tet_plane_inter_tris_batch

pp = np.array([0.5, 0., 0.])      # A point on the plane.
normal = np.array([1.,0.,0.])     # Vector normal to the plane.
tetrahedra = np.array([           # Tetrahedra to test for intersection.
    [[0.0696218 ,  0.51675287,  0.27622866],
     [0.77795185,  0.35896   ,  0.76889514],
     [0.93955095,  0.46952051,  0.50414596],
     [0.88091001,  0.0994195 ,  0.58077994]],
    [[0.22094245,  0.2076138 ,  0.97860943],
     [0.07205741,  0.20692974,  0.03908092],
     [0.42881122,  0.14083737,  0.73928276],
     [0.81604653,  0.96100499,  0.73425617]],
    [[0.36762177,  0.19788064,  0.32692846],
     [0.74746389,  0.26432934,  0.82284067],
     [0.59302365,  0.03474258,  0.48759837],
     [0.90889265,  0.50532202,  0.09855027]]
])

pint, area, center, tri_num, tri_ptr =\
        tet_plane_inter_tris_batch(pp, normal, tetrahedra)
```

The call returns an ndarray of intersection points (triangle vertex coords), an ndarray of triangle
areas, an ndarray of triangle center coordinates, the number of intersection triangles for each 
of the input tetrahedra, and the offsets of the data for each tetrahedron in the output arrays.
