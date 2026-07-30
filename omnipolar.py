import numpy as np


# *************************************************************************
# QCEP ITACA UPV
# Omnipolar Analysis
#
# Authors: Izan Segarra, Samuel Ruipérez-Campillo, Francisco Castells.
#          (Python port: Raúl Alós and EP Analytics Lab group)
#
# Any individual benefiting from any of this code must cite the work as:
# F. Castells, S. Ruipérez-Campillo, I. Segarra, R. Cervigón,
# R. Casado-Arroyo, J. Merino, J. Millet, Performance assessment of
# electrode configurations for the estimation of omnipolar electrograms
# from high density arrays, Computers in Biology and Medicine (2023).
#
# Description: Computes the orthogonal bipolar EGMs (bx, by) for every
#              clique of a RECTANGULAR high-density grid, using the
#              selected clique geometry.
# *************************************************************************
#
# COMPUTE_B_EGM computes the bipolar EGM on the x and y axes for each of the
# cliques of a rectangular grid using the specified method (Triangles, Mean
# Square or Cross).
#
#   b_egm = compute_b_egm(u_egm, n_rows, n_columns, mode)
#
#   Parameters:
#       u_egm (float ndarray): Unipolar electrograms, [n_electrodes x T],
#           ordered row by row (electrode (r,c) -> row (r-1)*n_columns + c).
#       n_rows (int): Number of rows of the catheter.
#       n_columns (int): Number of columns of the catheter.
#       mode (int): 1-4 Triangular, 5 Mean Square, 6 Cross.
#
#   Returns:
#       b_egm (float ndarray): Bipolar electrograms, [2*(n_rows-1)*(n_columns-1) x T],
#           interleaved per clique as [bx1; by1; bx2; by2; ...].
def compute_b_egm(u_egm, n_rows, n_columns, mode):
    u = np.asarray(u_egm, float)
    n_clique = (n_rows - 1) * (n_columns - 1)
    b_egm = np.zeros((2 * n_clique, u.shape[1]))
    Q45 = _rot2d(45)
    for row in range(1, n_rows):
        for col in range(1, n_columns):
            clique = (row - 1) * (n_columns - 1) + col          # 1-based
            e = (row - 1) * n_columns + col - 1                 # 0-based electrode
            nc = n_columns
            if mode == 1:
                bx = u[e + nc + 1] - u[e + nc]
                by = u[e] - u[e + nc]
            elif mode == 2:
                bx = u[e + 1] - u[e]
                by = u[e] - u[e + nc]
            elif mode == 3:
                bx = u[e + 1] - u[e]
                by = u[e + 1] - u[e + nc + 1]
            elif mode == 4:
                bx = u[e + nc + 1] - u[e + nc]
                by = u[e + 1] - u[e + nc + 1]
            elif mode == 5:
                bx = np.mean(np.vstack([u[e + 1] - u[e],
                                        u[e + nc + 1] - u[e + nc]]), 0)
                by = np.mean(np.vstack([u[e] - u[e + nc],
                                        u[e + 1] - u[e + nc + 1]]), 0)
            elif mode == 6:
                bx = u[e + 1] - u[e + nc]
                by = u[e] - u[e + nc + 1]
                bxy = Q45 @ np.vstack([bx, by])
                bx, by = bxy[0], bxy[1]
            else:
                raise ValueError('Introduced mode is not recognized.')
            c0 = clique - 1
            b_egm[2 * c0] = bx
            b_egm[2 * c0 + 1] = by
    return b_egm


# *************************************************************************
# EP Analytics Lab - ITACA UPV
# Omnipolar Analysis
#
# Authors: Raúl Alós and EP Analytics Lab group
#
# Any individual benefiting from any of this code must cite the work as:
# Alós R, Ramírez E, Peris A, Millet J, Castells F. Staggered arrangement of
# high-density multielectrodes for improved omnipolar sensing in cardiac
# electrophysiology: an in silico study. Computers in Biology and Medicine.
# 2026;214:111866. doi:10.1016/j.compbiomed.2026.111866
#
# Description: Computes the orthogonal bipolar EGMs (bx, by) for a clique of a
#              STAGGERED high-density grid (regular triangular, rhomboid or
#              hexagonal geometry).
# *************************************************************************
#
# COMPUTE_B_EGM_STAGGERED computes the orthogonal bipolar EGM (bx, by) of one
# staggered clique from its ordered node unipolar EGMs.
#
#   b_egm = compute_b_egm_staggered(u_clique, geometry, d)
#
#   Parameters:
#       u_clique (float ndarray): Ordered node unipolar EGMs of the clique,
#           [n_nodes x T] (3 for triangular, 4 for rhomboid, 7 for hexagonal).
#       geometry (str): 'triangular' | 'rhomboid' | 'hexagonal'.
#       d (float): Interelectrode distance (mm).
#
#   Returns:
#       b_egm (float ndarray): Orthogonal bipolar pair, [2 x T] = [bx; by].
def compute_b_egm_staggered(u_clique, geometry, d):
    N = np.asarray(u_clique, float)
    if geometry == 'triangular':
        bx = (N[1] - N[0]) / d
        by = ((N[2] - N[0]) + (N[2] - N[1])) / (np.sqrt(3) * d)
    elif geometry == 'rhomboid':
        bx = (N[1] - N[0]) / d
        by = (N[3] - N[2]) / (np.sqrt(3) * d)
    elif geometry == 'hexagonal':
        bx = ((N[3] - N[6]) / (3 * d)) + ((N[4] - N[5] + N[1] - N[0]) / (6 * d))
        by = ((N[5] - N[1]) + (N[4] - N[0])) / (2 * d * np.sqrt(3))
    else:
        raise ValueError('Introduced geometry is not recognized.')
    return np.vstack([bx, by])


# *************************************************************************
# EP Analytics Lab - ITACA UPV
# Omnipolar Analysis
#
# Authors: Raúl Alós and EP Analytics Lab group.
#
# Any individual benefiting from any of this code must cite the work as:
# Alós R, Ramírez E, Peris A, Millet J, Castells F. Staggered arrangement of
# high-density multielectrodes for improved omnipolar sensing in cardiac
# electrophysiology: an in silico study. Computers in Biology and Medicine.
# 2026;214:111866. doi:10.1016/j.compbiomed.2026.111866
#
# Description: General omnipolar computation. Given any set of orthogonal
#              bipolar pairs and the corresponding clique-centre (gold-standard)
#              unipolar EGMs, it rotates each pair onto its propagation
#              direction to obtain the omnipole, and returns the omnipolar EGM,
#              the Residual-to-Omnipolar Ratio (ROR) and the LAT error.
# *************************************************************************
#
# COMPUTE_O_EGM rotates the orthogonal bipolar pairs onto the propagation
# direction to obtain the omnipolar EGMs, and derives the ROR and LAT error.
#
#   o_egm, ror, lat_error, angle = compute_o_egm(b_egm, centroid_gt, fs)
#
#   Parameters:
#       b_egm (float ndarray): Orthogonal bipolar EGMs, [2*n_clique x T],
#           interleaved per clique as [bx1; by1; bx2; by2; ...].
#       centroid_gt (float ndarray): Clique-centre unipolar EGMs (gold standard
#           for local activation), [n_clique x T].
#       fs (float): Sampling rate of the EGMs (Hz). Default 1000*100.
#
#   Returns:
#       o_egm (float ndarray): Omnipolar EGMs, [2*n_clique x T], interleaved as
#           [main1; residual1; main2; residual2; ...].
#       ror (float ndarray): Residual-to-Omnipolar Ratio per clique, [n_clique].
#       lat_error (float ndarray): LAT error per clique (s), [n_clique].
#       angle (float ndarray): Estimated propagation angle per clique (rad,
#           in [0, 2*pi)), [n_clique].
def compute_o_egm(b_egm, centroid_gt, fs=1000 * 100):
    b_egm = np.asarray(b_egm, float)
    centroid_gt = np.atleast_2d(centroid_gt)
    n_clique = b_egm.shape[0] // 2
    T = b_egm.shape[1]
    o_egm = np.zeros((2 * n_clique, T))
    ror = np.zeros(n_clique)
    lat_error = np.zeros(n_clique)
    angle = np.zeros(n_clique)
    for c in range(n_clique):
        bx = b_egm[2 * c]
        by = b_egm[2 * c + 1]

        # propagation direction: angle of the largest-magnitude sample, [0, 2*pi)
        th = np.arctan2(by, bx)
        a = th[np.argmax(np.hypot(bx, by))]
        if a < 0:
            a += 2 * np.pi
        angle[c] = a

        # rotate onto the propagation direction -> omnipole (main) + residual
        R = np.array([[np.cos(-a), -np.sin(-a)], [np.sin(-a), np.cos(-a)]])
        oegm = R @ np.vstack([bx, by])
        o_egm[2 * c] = oegm[0]
        o_egm[2 * c + 1] = oegm[1]

        # metrics
        ror[c] = np.max(np.abs(oegm[1])) / np.max(np.abs(oegm[0]))
        t_gt = np.argmax(-np.gradient(centroid_gt[c]))
        t_omni = np.argmax(oegm[0])
        lat_error[c] = (t_gt - t_omni) / fs
    return o_egm, ror, lat_error, angle


def _rot2d(deg):
    r = np.radians(deg)
    return np.array([[np.cos(r), -np.sin(r)], [np.sin(r), np.cos(r)]])
