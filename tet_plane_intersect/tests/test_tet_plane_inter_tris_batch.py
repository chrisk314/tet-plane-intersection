import numpy as np

from tet_plane_intersect import tet_plane_inter_tris_batch


def test_plane_tet_inter_tris_batch():
    pp = np.array([0.5, 0., 0.])
    normal = np.array([1.,0.,0.])
    tets = np.array([
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
            tet_plane_inter_tris_batch(pp, normal, tets)

    assert len(area) == 3
    np.testing.assert_array_almost_equal(
        area, np.array([1.6705e-02, 6.2246e-02, 1.5645e-02])
    )

    assert center.shape == (3,3)
    np.testing.assert_array_almost_equal(
        center,
        np.array([
            [+5.000e-01, +4.03209e-01, +4.67449e-01],
            [+5.000e-01, +4.97728e-01, +6.80444e-01],
            [+5.000e-01, +1.98727e-01, +3.97374e-01]
        ])
    )

    np.testing.assert_array_equal(tri_num, np.array([1,1,1]))
    np.testing.assert_array_equal(tri_ptr, np.array([0,1,2,3]))
