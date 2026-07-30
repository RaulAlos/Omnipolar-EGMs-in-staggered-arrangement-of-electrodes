# Omnipolar EGMs in a staggered arrangement of electrodes

Simulation data and example accompanying the study:

> **Staggered arrangement of high-density multielectrodes for improved omnipolar sensing in cardiac electrophysiology: an in silico study.**
> Raúl Alós, Elisa Ramírez, Alfred Peris, José Millet, Francisco Castells.
> *Computers in Biology and Medicine* **214** (2026) 111866.
> DOI: [10.1016/j.compbiomed.2026.111866](https://doi.org/10.1016/j.compbiomed.2026.111866)

This repository provides the in silico simulation data (healthy and fibrotic
ventricular tissue) and a compact, self-contained script that reconstructs the
omnipolar electrogram (oEGM) for the five clique geometries evaluated in the
paper — Cross, Right Triangular, Regular Triangular, Rhomboid and Hexagonal —
on both a rectangular and a staggered electrode grid. It is an illustrative
example for reproducibility, not the full analysis pipeline of the study.

## What the example does

For a chosen tissue, anisotropy, lesion size, interelectrode distance and
catheter position/orientation, `main.py`:

1. loads the unified simulation data (surface node coordinates + electrograms);
2. shows the amplitude and local activation time (LAT) maps with both catheters
   overlaid (and the fibrotic lesion contours when applicable);
3. samples the extracellular field at several nearby random catheter positions
   ("freeze groups") and, for each clique geometry, reconstructs the local
   electric field, the omnipole and the residual;
4. reports the Residual-to-Omnipolar Ratio (ROR) and the LAT error, and plots
   the unipoles, bipoles, omnipole + residual and the bipolar loops.

## Repository contents

```
├── main.py                     Reproducibility example (entry point)
├── omnipolar.py                Bipolar / omnipolar EGM computation (ROR, LAT error)
├── functions.py                Geometry and field-sampling helpers
├── requirements.txt
├── LICENSE                     CC BY 4.0
└── CITATION.cff

Data (downloaded separately from Zenodo — see "Data" below):
├── Ventricular Healthy Data/   healthy_ar{03,05,07}_{coords,egms}.csv
└── Ventricular Fibrosis Data/  fibrotic_{size}_{R}_ar{aniso}_{coords,egms,
                                nucl_bound,fibro_bound_full}.csv
```

## Requirements

- Python 3.9+
- `numpy`, `scipy`, `matplotlib`

```
pip install -r requirements.txt
```

## Data

The electrogram data is archived on Zenodo:

**Download:** https://doi.org/10.5281/zenodo.21671892

Unzip it inside the repository folder so the two data folders sit next to
`main.py`:

```
Omnipolar-EGMs-in-staggered-arrangement-of-electrodes/
├── main.py
├── Ventricular Healthy Data/
└── Ventricular Fibrosis Data/
```

## Usage

Run the default case (healthy tissue, anisotropy 0.7, 5 mm spacing):

```
python main.py
```

Override any parameter from the command line, for example the fibrotic case:

```
python main.py --tissue fibrotic --model R2 --aniso 07
```

See all options with:

```
python main.py --help
```

| Argument      | Values                          | Meaning                                   |
|---------------|---------------------------------|-------------------------------------------|
| `--tissue`    | `healthy` \| `fibrotic`         | tissue type                               |
| `--aniso`     | `03` \| `05` \| `07`            | anisotropy ratio (0.3 / 0.5 / 0.7)        |
| `--model`     | `R1` \| `R15` \| `R2`           | fibrotic lesion size (2 / 3 / 4 cm)       |
| `--d`         | `2` \| `3` \| `4` \| `5`        | interelectrode distance (mm)              |
| `--Psi`       | degrees                         | catheter rotation                         |
| `--shift_x/y` | mm                              | catheter displacement                     |
| `--n_groups`  | integer                         | number of freeze groups                   |

Available fibrotic combinations: anisotropy `03/05/07` at `R2`, and lesion sizes
`R1/R15/R2` at anisotropy `07`.

## Data description

Each dataset is stored as plain-text CSV:

- `*_coords.csv` — surface node coordinates `[N x 3]` in micrometres (this is the
  surface mesh used in the study).
- `*_egms.csv` — extracellular potential at those nodes `[N x T]`.
- `*_nucl_bound.csv` / `*_fibro_bound_full.csv` — (fibrotic only) core and
  full-lesion boundary contours.

## Simulation methodology

Electrograms were generated with [openCARP](https://opencarp.org/) using a
pseudo-bidomain formulation. The repository shares the **superior (top) surface**
of the 3-D simulated tissue — the extracellular potential sampled on the upper
layer where the virtual electrodes are placed.

**Tissue and meshes.** Three 3-D slabs of healthy ventricular tissue
(Ten Tusscher–Panfilov ionic model), 8 × 8 × 0.5 cm, were immersed in a
conductive blood bath of 8 × 8 × 1 cm and meshed with tetrahedral elements at a
spatial resolution (node spacing) of 0.03 cm (300 µm). Node coordinates in the
`*_coords.csv` files are in micrometres on that surface.

**Anisotropy and conduction velocity.** The longitudinal conduction velocity
was fixed at 1.5 m/s, and the transverse/longitudinal conductivities were tuned
to anisotropy ratios ar = 0.3, 0.5 and 0.7, yielding three propagation patterns
of increasing wavefront planarity.

**Stimulation and recording.** A single 1 ms intracellular stimulus was applied
at t = 100 ms, eliciting one tissue activation. Extracellular potentials were
recorded at a temporal resolution of 1 ms (sampling frequency 1 kHz); the
healthy recordings span a single beat (1000 samples, ≈1 s). Unipolar EGMs were
obtained by spatially sampling the extracellular potential with 6 × 8 virtual
electrode arrays (rectangular and staggered).

**Fibrotic tissue.** Fibrosis was induced by ionic remodeling embedded at the
centre of the tissue as three spherical regions of 2, 3 and 4 cm diameter and
0.5 cm thickness (files tagged `R1`, `R15`, `R2`). The inner 75 % of each region
is the lesion core (homogeneous remodeling) and the outer 25 % the border zone
(healthy and remodeled elements in equal proportion); the `*_nucl_bound.csv` and
`*_fibro_bound_full.csv` files give the corresponding contours.

## Data availability

The code is hosted in this repository; the simulation data (electrograms) is
archived on Zenodo at https://doi.org/10.5281/zenodo.21671892 (CC BY 4.0).

## Citation

If you use this software or data, please cite the article above. A machine-readable
citation is provided in `CITATION.cff`.

## License

Code and data are released under the Creative Commons Attribution 4.0
International License (CC BY 4.0). See `LICENSE`.