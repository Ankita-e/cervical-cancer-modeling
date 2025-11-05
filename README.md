# Cervical Cancer Mortality Modelling (UK, 1993–2040)

This repository accompanies the research paper:

**Chakraborty, A. (2025).**  
*A Mathematical Modelling Framework for Cervical Cancer Mortality in the UK: Projections Towards 2040 Elimination.*



## Overview

This repository contains the MATLAB implementation of a mechanistic cervical cancer model, including:

- Parameter calibration using national incidence and mortality data (1993–2019)
- Optimal control analysis based on Pontryagin’s Maximum Principle
- Forward–Backward Sweep Method for computing optimal prevention and treatment controls
- Scenario forecasting to assess UK elimination feasibility by 2040
- Sensitivity analysis of intervention effectiveness

The model captures the dynamics of:
- **Incidence** \( I(t) \): declining through vaccination & screening
- **Case-Fatality Ratio** \( C(t) \): improving through treatment & early detection
- **Mortality** \( M(t) = I(t)*C(t) \)



## Repository Structure
cervical-cancer-modeling/
│
├── CODE/ # MATLAB scripts
│ ├── run_calibration_and_scenarios.m # Calibration, forecasting & sensitivity analysis
│ └── convergence.m # Forward–Backward Sweep optimal control solver
│
├── DATA/
│ └── cervical_cancer_tidy.xlsx # Cleaned UK national dataset (1993–2019)
│
├── FIGURE/ # Saved output figures (auto-generated)
│ ├── convergence_error.png # Convergence history of the FBSM method (used in paper)
│ ├── fit_plots.png # Observed versus model-fitted incidence and case-fatality ratio (CFR) (used in paper)
│ ├── states_control.png # Model states + optimal controls (used in paper)
│ ├── scenario_incidence.png # Projected incidence under intervention scenarios (used in paper)
│ └── sensitivity_gamma1.png # Effect of prevention strength on elimination year
│ └── scenario_cfr.png # (Generated, optional, not included in manuscript)
│
└── README.md


## Requirements

- MATLAB R2021a or later
- No external dependencies required



## How to Run

### Clone the repository:
```bash
git clone https://github.com/Ankita-e/cervical-cancer-modeling.git
cd cervical-cancer-modeling/CODE


## Citation
If using this repository, please cite:
Chakraborty, A. (2025). A Mathematical Modelling Framework for Cervical Cancer Mortality 
in the UK: Projections Towards 2040 Elimination.
