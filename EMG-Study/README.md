# Intra- and Inter-individual Spectral Pattern Variability of sEMG in Elbow Flexor Motor Tasks

## Overview

The pipeline performs:

1. Loading EMG recordings (`.mat`) from Noraxon Ultium, Qualisys/Delsys, or Neurosoft
2. Segmentation of continuous EMG into individual repetitions
3. Hann windowing and zero-padding to equal length
4. FFT-based spectral analysis
5. Spectral centroid computation
6. Condition comparison using Minkowski distance metrics
7. Generation of publication-ready figures (IEEE format)

Intermediate results are cached to avoid recomputation.

---

## Supported Muscles & Conditions

- **Muscles:**  
  - Brachioradialis (BR)  
  - Biceps brachii (BB)

- **Conditions:**  
  - Neutral grip  
  - Supinated grip

This repository contains data collected with Noraxon Ultium wireless hardware. 

---

## Structure
### Files
Go to project directory and run *main* script in console, just like below
```
main
```
This script takes raw EMG recordings, segments exercise repetitions, normalizes them, performs spectral and centroid analysis, compares conditions using multiple Minkowski distance formulations, and assembles publication-ready figures for an academic paper.

```
afterCalculations
```
This script is good to plot figures after all calculations are done. For free plotting when workspace is ready.
### Folders
All **.mat* files are in **data** folder

**processData** folder contains key functions

**additionalFiles** contains all rest files

---

#### Protocol how to collect data can by found here: [PL prototkół](protokół.md)
Project uses [fig Library](https://github.com/informacja/fig)
