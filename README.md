# Cervical Cancer Mortality Modelling (UK, 1993–2040)

This repository accompanies the research paper:

**Chakraborty, A. (2025).**  
*A Mathematical Modelling Framework for Cervical Cancer Mortality in the UK: Projections Towards 2040 Elimination.*

---

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
- **Mortality** \( M(t) = I(t) \cdot C(t) \)

---

## Repository Structure


