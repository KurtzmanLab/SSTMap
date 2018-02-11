import numpy as np

def rotate_check(matrix):

    if not (0.99 < np.linalg.det(matrix) < 1.01):

        raise Warning("Warning: Determinant of rotation matrix is %s. Should be close to +1.0." %np.linalg.det(matrix))


def do_rotation (crds, origin, rot_mat): 

    return (crds - origin).dot(rot_mat) + origin
