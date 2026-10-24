# MAST32001 Computational Statistics — Group O Final Project

Final assignment for MAST32001 (Computational Statistics), Group O:
Kaiwei Fu, Lingyu Wang, Sami Pönkä, Zhi Zhou.

Two tasks:

- **Task 1: Time series bootstrap and signal dimension.** Fixed-block and
  stationary bootstrap, the Ladle criterion, a simulation study, and an
  application to the course data (`timeseries.dat`).
- **Task 2: Bayesian mixture modelling with MCMC.** Two-component Gaussian
  mixture, a Gibbs sampler, label-switching analysis, and an application to
  the Old Faithful geyser data.

## Repository structure

```
zauziii-compstat-final-project/
├── R/                        # Analysis scripts
│   ├── project_setup.R       # Shared setup: project root, packages, theme, PNG export
│   ├── task1_functions.R     # Task 1 core functions
│   ├── task1_abc_analysis.R  # Task 1 Parts A-C: whitening, bootstrap, Ladle
│   ├── task1_de_analysis.R   # Task 1 Parts D-E: simulation study, course-data application
│   └── task2.R               # Task 2 Parts A-D: Gibbs sampler, simulations, Faithful
├── report/
│   └── refs.bib              # References
├── results/
│   ├── tables/               # All tables referenced by the report (CSV)
│   ├── data/                 # Intermediate outputs (checkpoints, criteria)
│   └── task2_sessionInfo.txt
└── archive/                  # Archived drafts and old outputs
```

Requires R packages `JADE` and `posterior`, and Quarto to render the PDF report.
