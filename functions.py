# *************************************************************************
# EP Analytics Lab - ITACA UPV
# Authors: Raúl Alós and EP Analytics Lab group
# Latest Update: 27/07/2026
#
# Helper functions for the omnipolar reconstruction example
# (Python port of angle_shift.m, applyRotationToMesh.m and gridOmnipolar.m,
#  plus movmean / field-sampling utilities used by main.py).
# *************************************************************************
import numpy as np
from scipy.spatial import Delaunay


def angle_shift(b_egm):
    """Dominant orientation of a 2xT bipolar pair, in [0, 2*pi)."""
    th = np.arctan2(b_egm[1, :], b_egm[0, :])
    r = np.hypot(b_egm[0, :], b_egm[1, :])
    a = th[np.argmax(r)]
    if a < 0:
        a += 2 * np.pi
    return a


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


def grid_omnipolar(cliques):
    """Reconstruct the omnipole and its metrics for the five clique geometries.

    cliques : dict with
        cliques['params']            -> [d, Psi, shift, aniso] (d used as step k)
        cliques[cath][geom]          -> dict with
            'node_pot'        (nC, nNodes, T)
            'centroid_pot'    (nC, T)
            'centroid_coords' (nC, 2)
      for (cath, geom) in optrell/SQ, optrell/TRI, stag/TRI, stag/RH, stag/HEX.

    Returns a dict with the same layout, each geometry holding
        'localEfield' (2*nC, T)  stacked [Exx; Eyy]
        'omnipolar'   (nC, 2, T) row 0 = omnipole, row 1 = residual
        'performance' (nC, 3)    [LAT_error (s), ROR, angle (rad)]
    """
    FS = 1000 * 100
    k = cliques['params'][0]
    efield = {
        'optrell,SQ':  lambda N: np.vstack([((N[2]-N[3])+(N[1]-N[0]))/(2*k),
                                            ((N[3]-N[0])+(N[2]-N[1]))/(2*k)]),
        'optrell,TRI': lambda N: np.vstack([(N[1]-N[2])/k, (N[2]-N[0])/k]),
        'stag,TRI':    lambda N: np.vstack([(N[1]-N[0])/k,
                                            ((N[2]-N[0])+(N[2]-N[1]))/(np.sqrt(3)*k)]),
        'stag,RH':     lambda N: np.vstack([(N[1]-N[0])/k, (N[3]-N[2])/(np.sqrt(3)*k)]),
        'stag,HEX':    lambda N: np.vstack([((N[3]-N[6])/(3*k)) + ((N[4]-N[5]+N[1]-N[0])/(6*k)),
                                            ((N[5]-N[1])+(N[4]-N[0]))/(2*k*np.sqrt(3))]),
    }
    geoms = [('optrell', 'SQ', 4), ('optrell', 'TRI', 3), ('stag', 'TRI', 3),
             ('stag', 'RH', 4), ('stag', 'HEX', 7)]
    results = {'optrell': {}, 'stag': {}, 'params': cliques['params']}
    for cath, geom, nN in geoms:
        clq = cliques[cath][geom]
        node_pot = np.asarray(clq['node_pot'], float)
        cen_pot = np.asarray(clq['centroid_pot'], float)
        cc = np.atleast_2d(clq['centroid_coords'])
        nC, T = cc.shape[0], cen_pot.shape[1]
        Exx, Eyy = np.zeros((nC, T)), np.zeros((nC, T))
        omni = np.zeros((nC, 2, T))
        ang, LATerr, ROR = np.zeros(nC), np.zeros(nC), np.zeros(nC)
        for i in range(nC):
            N = node_pot[i, :nN, :]
            uni = cen_pot[i, :]
            E = efield[f'{cath},{geom}'](N)
            Exx[i], Eyy[i] = E[0], E[1]
            phi = angle_shift(E)
            Rphi = np.array([[np.cos(-phi), -np.sin(-phi)],
                             [np.sin(-phi),  np.cos(-phi)]])
            oegm = Rphi @ E
            omni[i] = oegm
            t_uni = np.argmax(-np.gradient(uni))
            t_omni = np.argmax(oegm[0])
            LATerr[i] = (t_uni - t_omni) / FS
            ROR[i] = np.max(np.abs(oegm[1])) / np.max(np.abs(oegm[0]))
            ang[i] = phi
        results[cath][geom] = {'localEfield': np.vstack([Exx, Eyy]),
                               'omnipolar': omni,
                               'performance': np.column_stack([LATerr, ROR, ang])}
    return results
