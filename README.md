# Shared-State Inverse Synthesis of a Dual-Mode Conformal Metasurface

This repository contains the Mathematica 9 code, numerical data, and figures associated with the manuscript:

> M. R. Khodja and G. A. Alomar, “Shared-State Inverse Synthesis of Dual-Mode Conformal Metasurfaces for Low Scattering and Telemetry,” 2026.

The study considers a circular, 64-cell conformal impedance boundary with two electromagnetic operating modes. The low-scattering mode reduces monostatic scattering over several plane-wave incidences, whereas the telemetry mode forms a directive beam from an internal line source. Both modes use the same finite alphabet of eight passive impedance states, but each mode uses a different spatial codebook. Here, a *codebook* is the 64-entry assignment that specifies which impedance state is applied to each cell.

The main contribution is the formulation and solution of this shared-state, dual-mode inverse problem. Inverse-source null-space freedom is used as a design mechanism, but null-space exploitation by itself is not the claimed novelty.

## Scope of the repository

The calculations establish feasibility at the effective impedance-boundary level. They do not establish unit-cell, bias-network, switching-circuit, experimental, or platform-level realizability. The repository should therefore be interpreted as a reproducible computational study of the prescribed boundary model, not as an experimental validation package.

## Repository contents

- `RunDualModeRIS.nb`: launcher notebook for the primary numerical study.
- `DualModeRIS.m`: Mathematica 9 package containing the midpoint method-of-moments model, inverse synthesis, self-tests, parameter studies, data export, and plotting routines.
- `data/`: descriptive CSV files supporting Figs. 1–5 and Tables II–III of the manuscript.
- `figures/`: PDFs generated directly by `DualModeRIS.m`.
- `figures/Clean_Figs/`: a formatting notebook, local copies of the plotting data, and the reformatted manuscript figures. This notebook permits changes to fonts, colors, legends, and panel sizing without rerunning the numerical study.
- `FourierGalerkin_Mathematica9_Validation/`: an independent Fourier–Galerkin cross-formulation, its checks, data, figures, and supporting documentation.

The independent validation folder contains its own detailed files:

- `RunFourierGalerkinValidation.nb`: launcher for the cross-validation study.
- `FourierGalerkinValidation.m`: self-contained Fourier–Galerkin solver.
- `README_FIRST.md`: detailed execution instructions and expected numerical fingerprint.
- `FORMULATION_AND_CHECKS.md`: formulation, conventions, and definitions of the numerical checks.
- `IMPLEMENTATION_AUDIT.md`: implementation audit and interpretation of the cross-solver comparison.
- `fourier_galerkin_results/`: generated validation summary, CSV files, and PDF figures.

## Requirements

- Wolfram Mathematica 9. The packages were written to remain compatible with this version.
- A local directory for which Mathematica has write permission.
- A PDF viewer and spreadsheet program are useful for inspecting the exported results, but they are not required to run the calculations.

No Python installation, external electromagnetic solver, or nonstandard Mathematica package is required. Later Mathematica versions may also work, but the supplied notebooks were prepared and checked with Mathematica 9.

## Reproducing the primary numerical study

1. Download or clone the repository to a writable local directory.
2. Keep `RunDualModeRIS.nb` and `DualModeRIS.m` in the same directory.
3. Start Mathematica with a fresh kernel and open `RunDualModeRIS.nb`.
4. Select **Evaluation > Evaluate Notebook**, or evaluate the input cells from top to bottom with `Shift+Enter`.
5. Confirm that Mathematica prints `SelfTest PASSED.`
6. Allow `RunStudy[]` to finish. The calculation refreshes the CSV files in `data/`, exports Figs. 1–5 to `figures/`, and displays the figures in the launcher notebook.

The default production calculation uses 200 deterministic trials in the five-percent impedance-error study. A short workflow check may be performed with:

```Mathematica
RunStudy[MonteCarloTrials -> 20]
```

This reduced run is not a replacement for the reported production calculation. For a noninteractive export, use:

```Mathematica
RunStudy[DisplayFigures -> False]
```

To preserve the supplied output files while testing a rerun, direct the results to a separate folder:

```Mathematica
RunStudy[
  StudyDirectory -> FileNameJoin[{NotebookDirectory[], "rerun"}]
]
```

The reported codebooks are frozen. The reproducible default is therefore `RefineCodebooks -> False`.

### Primary-study numerical fingerprint

A successful default run should return values close to the following attached-run results. Small last-digit differences may occur across systems.

| Quantity | Expected value |
|---|---:|
| Validation worst protected-sector ratio | -11.575 dB |
| Validation mean protected-sector ratio | -13.089 dB |
| Telemetry two-dimensional directivity | 13.098 dB |
| Telemetry peak sidelobe level | -13.601 dB |
| Telemetry main-sector efficiency | 0.6491 |
| Shared-library mismatch, initial to final | 1.05846 to 0.009777 |

Do not use the production output if `SelfTest[]` fails.

## Reproducing the independent Fourier–Galerkin validation

1. Open `FourierGalerkin_Mathematica9_Validation/RunFourierGalerkinValidation.nb` in Mathematica 9.
2. Evaluate the notebook from top to bottom.
3. Confirm that Mathematica prints `FGSelfTest PASSED.` The launcher cancels the production calculation if this check fails.
4. Allow the four retained Fourier half-orders, $M=64,128,192,256$, to finish.
5. Inspect `fourier_galerkin_results/fg_validation_summary.txt`, followed by the CSV files and figures in the same results folder.

The default calculation uses 720 observation angles and 50-digit evaluation of sensitive Bessel–Hankel products before forming machine-precision system matrices. Depending on the processor, memory, and Mathematica installation, the full run may require several minutes to roughly half an hour.

At the final order, the supplied run gives a worst protected-sector ratio of approximately -11.36 dB, a mean ratio of -13.11 dB, a telemetry directivity of 13.02 dB, and a peak sidelobe level of -13.08 dB. The order-convergence, residual, and declared cross-solver checks all pass.

The Fourier–Galerkin solver does not load or call `DualModeRIS.m`. It uses a different expansion and testing procedure, analytic diagonalization of the circular single-layer operator, and exact Fourier coefficients of the piecewise-constant 64-cell profile. Exact equality with the midpoint-MoM results is neither imposed nor expected because the two formulations resolve the cell discontinuities and self-interaction differently.

## Replotting the manuscript figures without rerunning the solvers

Open:

```text
figures/Clean_Figs/Clean_Final_Figs_Based_on_DualModeRIS_and_FourierGalerkinValidation.nb
```

and evaluate it from top to bottom. The notebook uses the CSV files stored beside it and exports the reformatted PDFs to the same directory. This is the recommended route for changing plot fonts, colors, legends, image sizes, or panel spacing without repeating the numerical calculations.

The CSV files in the top-level `data/` directory are the direct scientific outputs of `DualModeRIS.m`. The corresponding files in `figures/Clean_Figs/` are plotting copies. In particular, the clean-figure copy of `fig1_configuration_codebooks.csv` contains additional spreadsheet helper columns used during figure preparation; its first ten columns contain the numerical data exported by the main package.

## Figure and table data map

| Manuscript item | Panel or content | Source CSV file |
|---|---|---|
| Fig. 1 | Shared impedance alphabet | `data/fig1_shared_impedance_states.csv` |
| Fig. 1(a) | Low-scattering codebook and cell coordinates | `data/fig1_configuration_codebooks.csv` |
| Fig. 1(b) | Telemetry codebook and cell coordinates | `data/fig1_configuration_codebooks.csv` |
| Fig. 2(a) | Normalized singular-value spectra | `data/fig2a_inverse_normalized_singular_spectra.csv` |
| Fig. 2(b) | Alternating-iteration history | `data/fig2b_inverse_iteration_history.csv` |
| Fig. 3(a) | Protected-sector ratio versus incidence angle | `data/fig3a_incidence_sweep.csv` |
| Fig. 3(b) | Normal-incidence scattering pattern | `data/fig3b_normal_incidence_scattering_pattern.csv` |
| Fig. 4(a) | Telemetry polar pattern | `data/fig4a_telemetry_polar_pattern.csv` |
| Fig. 4(b) | Telemetry control word | `data/fig4b_telemetry_control_word.csv` |
| Fig. 5(a) | Mesh-convergence study | `data/fig5a_mesh_convergence.csv` |
| Fig. 5(b) | Frozen-impedance frequency sweep | `data/fig5b_frozen_impedance_frequency_sweep.csv` |
| Fig. 6 | Fourier–Galerkin modal convergence | `FourierGalerkin_Mathematica9_Validation/fourier_galerkin_results/data/fg_validation_convergence.csv` |
| Table II | Forward-verified nominal performance | `data/table2_forward_verified_performance.csv` |
| Table III | Monte Carlo percentiles | `data/table3_monte_carlo_percentiles.csv` |
| Table III | Individual Monte Carlo trials | `data/table3_monte_carlo_trials.csv` |

All numerical CSV files use descriptive column headings. The far-field and convergence plots can therefore be regenerated independently of the solver notebooks.

## Model and numerical conventions

- Two-dimensional transverse-magnetic formulation with time dependence $\exp(-\mathrm{i}\omega t)$.
- Circular boundary of radius $a=2\lambda_0$, corresponding to a diameter of $4\lambda_0$.
- Sixty-four piecewise-constant controllable cells.
- One shared alphabet of eight passive normalized impedances, $Z/\eta$, and two mode-dependent 64-cell codebooks.
- Thirteen low-scattering incidences from $-30^\circ$ to $+30^\circ$ in $5^\circ$ increments.
- Internal telemetry line source at $(1.45\lambda_0,0)$.
- A 128-unknown synthesis model and an independent 256-unknown midpoint-MoM validation model, supplemented by mesh, frequency, and impedance-error studies.
- An independent Fourier–Galerkin calculation with $2M+1$ complex modes at each retained half-order $M$.

Additional definitions, normalization conventions, and acceptance criteria are given in the manuscript and in the documentation supplied with the Fourier–Galerkin package.

## Citation

If this repository contributes to published work, please cite the associated manuscript:

```bibtex
@unpublished{KhodjaAlomar2026SharedState,
  author = {Mohamed R. Khodja and Ghala A. Alomar},
  title  = {Shared-State Inverse Synthesis of Dual-Mode Conformal
            Metasurfaces for Low Scattering and Telemetry},
  year   = {2026},
  note   = {Manuscript submitted for publication}
}
```

The citation should be updated with the DOI or arXiv identifier when a permanent bibliographic record becomes available.

## Contact

For questions concerning the formulation or the numerical package, contact:

**Mohamed R. Khodja**  
Department of Electrical Engineering, College of Engineering  
Prince Mohammad Bin Fahd University, Khobar 31952, Saudi Arabia  
Email: `mkhodja@pmu.edu.sa`
