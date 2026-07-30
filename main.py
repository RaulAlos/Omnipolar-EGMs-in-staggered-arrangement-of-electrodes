# *************************************************************************
# EP Analytics Lab - ITACA UPV
#
# Authors: Raúl Alós and EP Analytics Lab group
# Latest Update: 27/07/2026
#
# Any individual benefiting from any of this code must cite the work as:
# Alós R, Ramírez E, Peris A, Millet J, Castells F.
# Staggered arrangement of high-density multielectrodes for improved
# omnipolar sensing in cardiac electrophysiology: an in silico study.
# Computers in Biology and Medicine. 2026;214:111866.
# doi:10.1016/j.compbiomed.2026.111866
#
# Description:
# Reproducibility example. For a chosen tissue, anisotropy, lesion
# size, interelectrode distance and catheter position/orientation, it loads
# the unified simulation data, overlays both catheters on the amplitude and
# LAT maps, and for the five clique geometries at several nearby random
# positions (freeze groups) computes the bipolar EGMs (compute_b_egm /
# compute_b_egm_staggered) and the omnipole, ROR and LAT error
# (compute_o_egm), then plots the unipoles, bipoles, omnipole+residual and
# bipolar loops.
#
# Configuration is set in the CONFIG block below (cfg.*).
#
# Requires: functions.py, omnipolar.py
# *************************************************************************
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from functions import apply_rotation_to_mesh, movmean, sample_field
from omnipolar import compute_b_egm, compute_b_egm_staggered, compute_o_egm

# ---- CONFIGURATION ------------------------------------------------------
# Defaults below; override from the terminal, e.g.:
#   python main.py --tissue fibrotic --model R2 --aniso 07
p = argparse.ArgumentParser(description='Minimal omnipolar reproducibility example.')
p.add_argument('--tissue', choices=['healthy', 'fibrotic'], default='healthy')
p.add_argument('--aniso',  choices=['03', '05', '07'], default='07')
p.add_argument('--model',  choices=['R1', 'R15', 'R2'], default='R2')  # fibrotic lesion size
p.add_argument('--d',        type=int,   default=3)   # interelectrode distance (mm)
p.add_argument('--Psi',      type=float, default=0)  # catheter rotation (deg)
p.add_argument('--shift_x',  type=float, default=0)   # displacement X (mm)
p.add_argument('--shift_y',  type=float, default=0)   # displacement Y (mm)
p.add_argument('--n_groups', type=int,   default=5)   # number of freeze groups
p.add_argument('--fs_native', type=int,  default=1000)  # egms sampling rate (Hz)
cfg = vars(p.parse_args())
nRows, nCols = 6, 8

# ---- LOAD DATA (CSV) ----------------------------------------------------
if cfg['tissue'] == 'healthy':
    folder = 'Ventricular Healthy Data'
    base = f"healthy_ar{cfg['aniso']}"
    fibrosis = None
else:
    folder = 'Ventricular Fibrosis Data'
    szmap = {'R1': '2cm', 'R15': '3cm', 'R2': '4cm'}
    base = f"fibrotic_{szmap[cfg['model']]}_{cfg['model']}_ar{cfg['aniso']}"
    fibrosis = dict(
        nucl_bound=np.loadtxt(os.path.join(folder, f"{base}_nucl_bound.csv"), delimiter=','),
        fibro_bound_full=np.loadtxt(os.path.join(folder, f"{base}_fibro_bound_full.csv"), delimiter=','),
    )
coords = np.loadtxt(os.path.join(folder, f"{base}_coords.csv"), delimiter=',')
egms = np.loadtxt(os.path.join(folder, f"{base}_egms.csv"), delimiter=',')
cx, cy = coords[:, 0].mean(), coords[:, 1].mean()

# ---- ELECTRODE GEOMETRY -------------------------------------------------
d = cfg['d']
XO, YO = np.meshgrid(np.arange(nCols), np.arange(nRows))
XO, YO = XO * d * 1000.0, YO * d * 1000.0
XO = XO + cx - (XO.max() + XO.min()) / 2
YO = YO + cy - (YO.max() + YO.min()) / 2
XS, YS = np.meshgrid(np.arange(nCols), np.arange(nRows))
XS, YS = XS * d * 1000.0, YS * d * 1000.0 * np.sqrt(3) / 2
XS = XS + cx - (XS.max() + XS.min()) / 2
XS[1::2, :] = XS[1::2, :] + d * 1000 / 2
YS = YS + cy - (YS.max() + YS.min()) / 2
psi = np.radians(cfg['Psi'])
Rm = np.array([[np.cos(psi), -np.sin(psi)], [np.sin(psi), np.cos(psi)]])
XO, YO = apply_rotation_to_mesh(XO, YO, Rm)
XS, YS = apply_rotation_to_mesh(XS, YS, Rm)
XO, YO = XO + cfg['shift_x'] * 1000, YO + cfg['shift_y'] * 1000
XS, YS = XS + cfg['shift_x'] * 1000, YS + cfg['shift_y'] * 1000
xo, yo = XO.flatten('F'), YO.flatten('F')     # column-major = paper node order
xs, ys = XS.flatten('F'), YS.flatten('F')

# name, grid, method, node indices (0-based)
#   rect: node order = 2x2 cell [top-left, top-right, bottom-left, bottom-right]
#         method = compute_b_egm mode (5 mean-square, 6 cross, 1-4 triangular)
#   stag: method = geometry name for compute_b_egm_staggered
specs = [
    ('Cross',              'rect', 5,            np.array([1, 7, 2, 8]) - 1),
    ('Right Triangular',   'rect', 2,            np.array([1, 7, 2, 8]) - 1),
    ('Regular Triangular', 'stag', 'triangular', np.array([1, 7, 2]) - 1),
    ('Rhomboid',           'stag', 'rhomboid',   np.array([2, 8, 7, 9]) - 1),
    ('Hexagonal',          'stag', 'hexagonal',  np.array([7, 13, 8, 14, 15, 9, 2]) - 1),
]
nGeo = len(specs)

# freeze groups: nearby random offsets, same for both grids
radius_mm = cfg['d']
ang_g = 2 * np.pi * np.random.rand(cfg['n_groups'])
rad_g = radius_mm * np.sqrt(np.random.rand(cfg['n_groups']))
off = np.column_stack([rad_g * np.cos(ang_g), rad_g * np.sin(ang_g)]) * 1000

mk = lambda x, idx: np.append(x[idx], x[idx].mean())   # node coords + centroid
grids = {'rect': (xo, yo), 'stag': (xs, ys)}
inst = []
for g in range(cfg['n_groups']):
    for name, grid, method, idx in specs:
        gx, gy = grids[grid]
        inst.append(dict(group=g, name=name, grid=grid, method=method,
                         x=mk(gx, idx) + off[g, 0], y=mk(gy, idx) + off[g, 1]))

# ---- MAPS ---------------------------------------------------------------
amp = egms.max(1) - egms.min(1)
LAT = np.array([np.argmin(np.gradient(egms[i])) for i in range(egms.shape[0])])

fig, axs = plt.subplots(1, 2, figsize=(12, 5), num='Tissue maps')
for ax, (val, ttl) in zip(axs, [(amp, f"Amplitude map ({cfg['tissue']}, ar={cfg['aniso']})"),
                                 (LAT, 'LAT map (sample index)')]):
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=val, s=14)
    if fibrosis is not None:
        ax.plot(fibrosis['nucl_bound'][:, 0], fibrosis['nucl_bound'][:, 1], 'k-', lw=1.5)
        ax.plot(fibrosis['fibro_bound_full'][:, 0], fibrosis['fibro_bound_full'][:, 1], 'k--', lw=1.2)
    ax.plot(xo, yo, 'ko', mfc='k', ms=3, label='Rectangular cath.')
    ax.plot(xs, ys, 'o', color=(0.85, 0, 0.2), mfc=(0.85, 0, 0.2), ms=3, label='Staggered cath.')
    ax.set_aspect('equal'); ax.grid(True); ax.legend(loc='best')
    ax.set_title(ttl); ax.set_xlabel('X (µm)'); ax.set_ylabel('Y (µm)')
    fig.colorbar(sc, ax=ax)
fig.tight_layout()

# ---- SAMPLE FIELD -------------------------------------------------------
query, rng_c = [], []
for k in inst:
    n0 = len(query)
    query.extend(np.column_stack([k['x'], k['y']]).tolist())
    rng_c.append((n0, len(query)))
query = np.array(query)

P = sample_field(coords, egms, query)
T = egms.shape[1]
xo_old, xo_new = np.linspace(0, 1, T), np.linspace(0, 1, T * 100)
P = np.array([np.interp(xo_new, xo_old, P[i]) for i in range(P.shape[0])])
P = movmean(P, 200, axis=1)
for j, k in enumerate(inst):
    a, b = rng_c[j]
    k['pot'] = P[a:b, :]

# ---- BIPOLES -> OMNIPOLE (per freeze group) -----------------------------
for g in range(cfg['n_groups']):
    idxg = [j for j, k in enumerate(inst) if k['group'] == g]

    b_stack, gt_stack = [], []
    for j in idxg:
        k = inst[j]
        nodes = k['pot'][:-1, :]        # node unipolar EGMs
        gt = k['pot'][-1, :]            # clique-centre unipole (LAT gold standard)
        if k['grid'] == 'rect':
            b = compute_b_egm(nodes, 2, 2, k['method'])         # [2 x T]
        else:
            b = compute_b_egm_staggered(nodes, k['method'], d)  # [2 x T]
        k['Ef'] = b
        b_stack.append(b)
        gt_stack.append(gt)

    b_egm = np.vstack(b_stack)          # [2*nGeo x T]
    centroid_gt = np.vstack(gt_stack)   # [nGeo x T]
    o_egm, ror, lat, ang = compute_o_egm(b_egm, centroid_gt, fs=cfg['fs_native'] * 100)

    for i, j in enumerate(idxg):
        k = inst[j]
        k['oegm'] = o_egm[2 * i:2 * i + 2, :]   # [main; residual]
        k['ROR'] = ror[i]
        k['eLAT'] = lat[i] * 1000               # ms

# ---- SIGNALS ------------------------------
fs = cfg['fs_native'] * 100
tvec = np.arange(inst[0]['pot'].shape[1]) / fs * 1000
kref = next(k for k in inst if k['group'] == 0 and k['name'] == 'Cross')
iact = np.argmin(np.gradient(kref['pot'][-1, :]))
xwin = (max(tvec[0], tvec[iact] - 100), min(tvec[-1], tvec[iact] + 200))

for g in range(cfg['n_groups']):
    byname = {k['name']: k for k in inst if k['group'] == g}
    fig, ax = plt.subplots(nGeo, 4, figsize=(15, 9), num=f'Freeze group {g + 1}')
    for s, (name, grid, method, idx) in enumerate(specs):
        k = byname[name]
        last = (s == nGeo - 1)

        ax[s, 0].plot(tvec, k['pot'][:-1, :].T, lw=0.5)
        ax[s, 0].plot(tvec, k['pot'][-1, :], 'k--', lw=1.2, label='centroid')
        ax[s, 0].set_xlim(xwin); ax[s, 0].grid(True)
        ax[s, 0].set_ylabel(f"{name}\nROR = {k['ROR']:.3f}\neLAT = {k['eLAT']:.2f} ms",
                            fontweight='bold')
        if s == 0: ax[s, 0].set_title('Unipoles'); ax[s, 0].legend(loc='best', fontsize=7)
        if last: ax[s, 0].set_xlabel('Time (ms)')

        ax[s, 1].plot(tvec, k['Ef'][0], 'k', label='E_xx')
        ax[s, 1].plot(tvec, k['Ef'][1], 'k--', label='E_yy')
        ax[s, 1].set_xlim(xwin); ax[s, 1].grid(True)
        if s == 0: ax[s, 1].set_title('Bipoles'); ax[s, 1].legend(loc='best', fontsize=7)
        if last: ax[s, 1].set_xlabel('Time (ms)')

        ax[s, 2].plot(tvec, k['oegm'][0], 'k', label='omnipole')
        ax[s, 2].plot(tvec, k['oegm'][1], color=(.7, .7, .7), label='residual')
        ax[s, 2].set_xlim(xwin); ax[s, 2].grid(True)
        if s == 0: ax[s, 2].set_title('Omnipole + Residual'); ax[s, 2].legend(loc='best', fontsize=7)
        if last: ax[s, 2].set_xlabel('Time (ms)')

        ax[s, 3].plot(k['Ef'][0], k['Ef'][1], 'k')
        m = 1.05 * np.max(np.abs(np.concatenate([k['Ef'][0], k['Ef'][1]])))
        ax[s, 3].set_xlim(-m, m); ax[s, 3].set_ylim(-m, m)
        ax[s, 3].set_aspect('equal', adjustable='box')
        ax[s, 3].grid(True); ax[s, 3].set_ylabel('E_yy')
        if s == 0: ax[s, 3].set_title('Loop (E_xx vs E_yy)')
        if last: ax[s, 3].set_xlabel('E_xx')
    fig.tight_layout()

plt.show()
