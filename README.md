# SSA Noise Toy Model — CMS Phase-2 PS Modules

A Python-based analytical/statistical toy model for simulating correlated 
noise occupancy in Short Strip ASIC (SSA) chips of CMS Phase-2 Outer Tracker 
Pixel-Strip (PS) modules, designed to be interfaced with CMS official 
simulation at the digitization step.

## 📖 Documentation & Step-by-Step Instructions
 
For a full description of the workflow — including how to extract thresholds,
compute correlation matrices, estimate f_corr, run the model, and validate
against real data — please refer to the
**[Wiki](https://github.com/aman0404/PSmoduleNoiseModeling/wiki/SSA-Noise-Toy-Model-%E2%80%94-Wiki)**.

# Motivation

In the CMS Phase-2 simulation chain, the digitization step converts 
simulated energy deposits into detector signals and adds noise. Standard 
digitizers typically treat channel noise as uncorrelated. However, SSA chips 
in PS module hybrids exhibit significant correlated noise structures — both 
within a chip (common-mode across subchip regions) and across chips on the 
same hybrid.

This toy model captures these correlations analytically and is intended to 
be stitched into the CMS digitization framework, providing a more realistic 
noise background for tracking and trigger performance studies.

## Overview

Each PS module hybrid carries 8 SSA chips, each reading out 120 channels. 
This model simulates noise hits across all chips and channels by combining:

- **Per-channel threshold distributions** — loaded from real or derived 
  threshold calibration data (e.g. from occupancy-based threshold extraction).
- **Correlated common-mode noise** — modeled via a 16×16 subchip correlation 
  matrix. Each SSA chip is split into two subchips of 60 channels each, 
  yielding 16 subchip units per hybrid. Correlations between all pairs are 
  encoded in a symmetric matrix and decomposed via Cholesky factorization 
  to generate physically consistent correlated noise offsets per event.
- **Per-chip correlation fraction (f_corr)** — a chip-level parameter that 
  controls the mixing between uncorrelated channel noise and the shared 
  common-mode offset for each subchip region.

The noise for each channel is modeled as:

    noise = sqrt(1 - f_corr) * uncorrelated_noise + sqrt(f_corr) * subchip_offset

A hit is registered when the noise exceeds the channel's threshold.

## Features

- Simulation of noise hits over configurable numbers of events (default: 50,000)
- 16×16 inter-subchip correlation matrix support (Cholesky decomposition)
- Per-chip f_corr parameter loaded from external calibration files
- Per-channel threshold loading from calibration data
- Per-chip hit count histograms (PNG output, log scale)
- ROOT output: per-chip `TH1F` histograms and a `TTree` with per-event 
  hit counts for all chips
- Summary statistics (mean and sigma of hits per chip) printed to stdout

## Input Files

| File | Description |
|------|-------------|
| `Hybrid1_fcorr.txt` | Per-chip correlation fraction f_corr (dict: chip → float) |
| `ThresholdsFromOcc_Hyb1_2sigma.txt` | Per-channel thresholds (dict: chip → list of 120 values) |
| `SubchipCorrelationCoefficients_16x16_Hyb1.txt` | 16×16 subchip correlation matrix (dict: "i-j" → rho) |

## Output

- `NoiseHits_Chip{N}.png` — per-chip noise hit distributions
- `ToyHits_Hyb1_2sigma_subChip_f_corr.root` — ROOT file containing:
  - `TH1F` histograms per chip
  - `toyTree` TTree with per-event hit counts for all 8 chips

## Requirements

- Python 3.x
- NumPy
- Matplotlib
- ROOT / PyROOT

## Usage

Edit the file paths in `__main__` to point to your input files, then run:

    python new_toymodel_mixed_subChips_fcorr.py

## Context

This model is developed in the context of the CMS Phase-2 Outer Tracker 
upgrade. It is intended for threshold optimisation studies and for 
understanding the impact of correlated noise on the noise occupancy of 
PS module hybrids.

## What the Model Adds to Digitization

- Replaces the uncorrelated noise assumption in the SSA 
  digitizer with a physically motivated correlated noise model
- Allows threshold-dependent noise occupancy to be injected at the 
  digitization step using real calibration data
- Enables studies of how correlated noise affects strip hit rates, 
  cluster sizes, and ultimately tracking efficiency and fake rate in 
  the Phase-2 Outer Tracker
