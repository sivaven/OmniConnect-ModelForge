# Cross-Species Connectomics

This project is part of the **OmniConnect-ModelForge** umbrella repository.  
It implements a framework for constructing **directed connectomes** across mouse, marmoset, rhesus macaque, and human, integrating tracer-derived projection polarity with species-specific diffusion MRI (dMRI) tractography.  

---

## Project Structure

```
projects/cross_species_connectomics/
├── README.md                  # Project-specific documentation (this file)
├── requirements.txt           # Project-specific dependencies
│
├── cross_modal_validation.py  # Figure 2c: Validation analyses across tracer vs dMRI modalities
├── cross_species_nets.py      # Directed connectome construction, path efficiency, network utilities
├── stat_utils.py              # Bootstrap, non-parametric tests, statistical utilities
├── vis_utils.py               # Plotting utilities, visualization helpers
├── voxel_count_stats.py       # Gray matter voxel statistics, region counts
│
├── fig3a.py                   # Figure 3a generation script
├── fig3b-4a.py                # Figures 3b and 4a generation script
├── fig4c.py                   # Figure 4c generation script
├── fig5a.py                   # Figure 5a generation script
├── figs_aba_scale.py          # Figure 2d and 3d generation script
│
├── dsi_aba_mouse.sh           # Tractography: seed-only + seed→matched projection (all injections)
├── dsi_aba_mouse_perm.sh      # Tractography: seed→randomized projection (permutation control)
│
└── data/                      # Expected location for injection/projection density NIfTIs, connectomes, etc.
```

---

## Requirements

- **Python** ≥ 3.8  
- **Dependencies** (install via `pip install -r requirements.txt`):
  - `numpy`
  - `pandas`
  - `networkx`
  - `matplotlib`
  - `seaborn`
  - `scipy`
  - `statsmodels`
- **DSI Studio** (latest build) — required for tractography (shell scripts).

---

## Usage

### 1. Tractography

Run tractography for all 1,199 Allen mouse injections:

- **Matched (real) connectome**
  ```bash
  bash dsi_aba_mouse.sh
  ```
  Generates:
  - `<ID>_seed.tt.gz` (seed-only streamlines)
  - `<ID>_seed__<ID>_end.tt.gz` (seed → matched projection)

- **Permutation (control) connectome**
  ```bash
  bash dsi_aba_mouse_perm.sh 
  ```
  Generates:
  - `<ID>_seed__<randomID>_end.tt.gz` (seed → randomized projection)  

Outputs are saved in `trk_outputs/`.

---

### 2. Python Analysis & Figures

Each `fig*.py` file generates the corresponding manuscript figure. Example:

```bash
python fig3a.py
```

Figures and tables will be written to the working directory or specified output paths.

Utility modules:
- `cross_species_nets.py`: core connectome building and path efficiency metrics.  
- `stat_utils.py`: bootstrapping, Kruskal–Wallis, FDR correction, etc.  
- `vis_utils.py`: plotting helpers for connectome evolution and comparisons.  
- `voxel_count_stats.py`: voxel-based region metrics.  

> All scripts are wrapped with `if __name__ == "__main__":`, so modules can be imported without side effects.

---

### 3. Data Organization

Expected inputs:
- `injection_densities/` — Allen Mouse Brain Atlas injection density maps (NIfTI).  
- `projection_densities/` — Allen Mouse Brain Atlas projection density maps (NIfTI).  
- Connectome matrices, label lists, and supplementary CSVs in `data/`.  

Ensure paths inside the scripts (`_dir` variables) are updated to match your environment.

---

## Reproducibility

- All statistical tests use functions in `stat_utils.py` with explicit random number generators (`numpy.random.default_rng`).  
- Figures are generated directly from raw data → connectomes → plots.

---

## Citation

If you use this project, please cite:

> **Cross-species brain circuitry from diffusion MRI tractography and mouse viral tracing.**  
> Siva Venkadesh, Wen-Jieh Linn, Yuhe Tian, G Allan Johnson, Fang-Cheng Yeh.  
> *bioRxiv* (2025). doi: [10.1101/2025.09.07.674762](https://doi.org/10.1101/2025.09.07.674762)

---
