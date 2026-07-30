# *************************************************************************
# EP Analytics Lab - ITACA UPV
# Authors: Raúl Alós and EP Analytics Lab group
# Latest Update: 27/07/2026
#
# Geometry and sampling helpers used by main.py: mesh rotation, moving
# average and field interpolation.
# *************************************************************************
import numpy as np
from scipy.spatial import Delaunay


def apply_rotation_to_mesh(x_mat, y_mat, R):
    """Rotate a 2-D mesh about its own centroid."""
    cx, cy = x_mat.mean(), y_mat.mean()
    pts = np.vstack([x_mat.ravel() - cx, y_mat.ravel() - cy])
    rot = R @ pts
    return rot[0].reshape(x_mat.shape) + cx, rot[1].reshape(y_mat.shape) + cy


def movmean(x, k, axis=-1):
    """Centred moving average with shrinking edges (like MATLAB movmean)."""
    x = np.moveaxis(np.asarray(x, float), axis, -1)
    n = x.shape[-1]
    hb, hf = k // 2, k - k // 2 - 1
    cs = np.concatenate([np.zeros(x.shape[:-1] + (1,)), np.cumsum(x, -1)], -1)
    i = np.arange(n)
    lo, hi = np.clip(i - hb, 0, n), np.clip(i + hf + 1, 0, n)
    out = (cs[..., hi] - cs[..., lo]) / (hi - lo)
    return np.moveaxis(out, -1, axis)


def sample_field(coords, egms, query_xy):
    """Linear (barycentric) interpolation of egms at query points.
    NaN outside the convex hull (equivalent to scatteredInterpolant 'none').
    Note: MATLAB used 'natural'; 'linear' is a close approximation."""
    tri = Delaunay(coords[:, :2])
    simplex = tri.find_simplex(query_xy)
    P = np.full((query_xy.shape[0], egms.shape[1]), np.nan)
    inside = simplex >= 0
    if np.any(inside):
        Tr = tri.transform[simplex[inside]]
        dxy = query_xy[inside] - Tr[:, 2, :]
        b = np.einsum('ijk,ik->ij', Tr[:, :2, :], dxy)
        bary = np.hstack([b, 1 - b.sum(1, keepdims=True)])
        verts = tri.simplices[simplex[inside]]
        Pin = np.zeros((int(inside.sum()), egms.shape[1]))
        for j in range(3):
            Pin += bary[:, [j]] * egms[verts[:, j], :]
        P[inside] = Pin
    return P
