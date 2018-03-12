
"""Calculates intersections between tetrahedra and planes in 3D."""

import ctypes
import os

import numpy as np
from numpy import ctypeslib as ctl

def wrapped_ndptr(*args, **kwargs):
    base = ctl.ndpointer(*args, **kwargs)

    def from_param(cls, obj):
        if obj is None:
            return obj
        return base.from_param(obj)

    return type(
        base.__name__, (base,),
        {'from_param': classmethod(from_param)}
    )

_libwrap_tetinter = ctl.load_library(
    'libwrap_tetinter', os.path.join(os.path.dirname(__file__), 'src')
)

Float64ArrayType = wrapped_ndptr(dtype=np.float64, flags='aligned, c_contiguous')
Int32ArrayType = wrapped_ndptr(dtype=np.int32, flags='aligned, c_contiguous')

_libwrap_tetinter.plane_normal_tetrahedron_intersect.argtypes = [
    Float64ArrayType,
    Float64ArrayType,
    Float64ArrayType,
    Float64ArrayType,
    ctypes.POINTER(ctypes.c_int),
]

_libwrap_tetinter.batch_intersect_tris.argtypes = [
    Float64ArrayType,
    Float64ArrayType,
    Float64ArrayType,
    ctypes.c_int,
    Float64ArrayType,
    Float64ArrayType,
    Float64ArrayType,
    Int32ArrayType,
    Int32ArrayType,
]

_libwrap_tetinter.quad_area_3d.argtypes = [
    Float64ArrayType,
]
_libwrap_tetinter.quad_area_3d.restype = ctypes.c_double

_libwrap_tetinter.tri_area_3d.argtypes = [
    Float64ArrayType,
]
_libwrap_tetinter.tri_area_3d.restype = ctypes.c_double


def tet_plane_inter_points(pp, normal, tet):
    """Calculate intersection points for a single tetrahedron.
    :param pp: ndarray specifying a single point in R^3 on the intersection plane.
    :param normal: ndarray specifying vector in R^3 normal to the intersection plane.
    :param tet: ndarray specifying the 4 corners of a tetrahedron in R^3.
    :returns: ndarray specifying 0, 1, 2, 3, or 4 intersection points.
    """
    global _libwrap_tetinter

    pint = np.empty((4,3), dtype=np.float64)
    int_num = ctypes.c_int()

    _libwrap_tetinter.plane_normal_tetrahedron_intersect(
        np.ascontiguousarray(pp, np.float64),
        np.ascontiguousarray(normal, np.float64),
        np.ascontiguousarray(tet, np.float64),
        np.ascontiguousarray(pint, np.float64),
        ctypes.byref(int_num)
    )

    return pint[:int_num.value]


def tet_plane_inter_tris_batch(pp, normal, tets):
    """Calculate intersection triangles for a set of tetrahedra.
    :param pp: ndarray specifying a single point in R^3 on the intersection plane.
    :param normal: ndarray specifying vector in R^3 normal to the intersection plane.
    :param tets: ndarray specifying the 4 corners of n tetrahedra in R^3.
    :returns: tuple containing ndarrays specifying 0, 1, 2, 3, or 4 intersection
        points, areas, and center coordinates of each intersection triangle, as well
        as the number of intersection triangles and the offset of the data for each
        tetrahedron in the input.
    """
    global _libwrap_tetinter

    # Preallocate arrays
    n_tets = len(tets)
    pint = np.empty((6*n_tets,3), dtype=np.float64)
    int_area = np.empty(2*n_tets, dtype=np.float64)
    int_center = np.empty((2*n_tets,3), dtype=np.float64)
    tri_num = np.empty(n_tets, dtype=np.int32)
    tri_ptr = np.empty(n_tets+1, dtype=np.int32)

    # Call C level function
    _libwrap_tetinter.batch_intersect_tris(
        np.ascontiguousarray(pp, np.float64),
        np.ascontiguousarray(normal, np.float64),
        np.ascontiguousarray(tets, np.float64),
        n_tets,
        np.ascontiguousarray(pint, np.float64),
        np.ascontiguousarray(int_area, np.float64),
        np.ascontiguousarray(int_center, np.float64),
        np.ascontiguousarray(tri_num, np.int32),
        np.ascontiguousarray(tri_ptr, np.int32),
    )

    # Slice arrays to fit data
    n_tris = tri_ptr[-1]
    pint = pint[:3*n_tris].reshape((n_tris,3,3))
    int_area = int_area[:n_tris]
    int_center = int_center[:n_tris]

    return pint, int_area, int_center, tri_num, tri_ptr


def inter_poly_area(pint):
    """Return area for intersection polygon.
    :param pint: ndarray containing intersection points.
    :return: float area of polygon.
    """
    global _libwrap_tetinter

    if len(pint) == 4:
        area = _libwrap_tetinter.quad_area_3d(
            np.ascontiguousarray(pint, np.float64)
        )
    elif len(pint == 3):
        area = _libwrap_tetinter.tri_area_3d(
            np.ascontiguousarray(pint, np.float64)
        )
    else:
        area = 0.
    return area


def inter_poly_center(pint):
    """Return center of mass of intersection points.
    :param pint: ndarray containing intersection points.
    :return: ndarray containing COM of intersection points in R^3.
    """
    return np.mean(pint, axis=0)


class PlaneTetInter(object):

    """Represents the intersection of a plane and a tetrahedron."""

    def __init__(self, pp, normal, tet):
        self.points = tet_plane_inter_points(pp, normal, tet)
        self.n_points = len(self.points)
        self.area = inter_poly_area(self.points)
        if self.n_points > 0:
            self.center = inter_poly_center(self.points)
        else:
            self.center = None
