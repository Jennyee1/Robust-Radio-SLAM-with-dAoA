# Robust-Radio-SLAM-with-dAoA

This repository contains the official MATLAB implementation for the paper:

**"ROBUST RADIO SLAM BY FUSING TEMPORAL-SPATIAL DIFFERENTIAL AOA INFORMATION COOPERATIVELY"**

**Status:** Submitted to **ICASSP 2026**
> *Update (2025-12-21): The codebase has been updated to enhance the reproducibility of the new version. Detailed usage instructions have been added.*
---

## Introduction
This project implements a robust Radio-SLAM framework that cooperatively fuses **Temporal Differential AoA (T-dAoA)** and **Spatial Differential AoA (S-dAoA)** information. 

### Key Features
* **Temporal & Spatial Fusion:** Jointly utilizes time-domain and space-domain differential angle information.
* **Robust Data Association:** BP-based message passing for accurate feature identification in cluttered environments.
* **Detailed Analysis:** Includes ablation studies and complexity analysis (see Supplementary Material).

---

## Repository Structure

* **`main_example.m`**: The main entry script. Run this file to reproduce the simulation results and visualize the SLAM process.
* **`BPbasedSLAM_dAoA_v2.m`**: The core function implementing the proposed SLAM framework.
* **`functions/`**: Contains all necessary sub-functions for data generation, particle filtering, and mathematical calculations, etc.
* **`Get_ablation_res.m`** & **`Res_ablation.mat`**: Scripts and data files related to the ablation studies presented in the paper.
* **`SupplementaryMaterial.pdf`**: Provides detailed supplementary content, including:
    * Detailed Ablation Study Results
    * Computational Complexity Analysis
    * List of Symbols and Notations

---

## Requirements

* **MATLAB R2024b** (or later recommended)
* *Note: Ensure standard toolboxes

---

## Run the Simulation
1. Open MATLAB and navigate to the repository folder.
2. Add the functions folder to your MATLAB path (or let the script handle it).
3. Run the main example script: main_example.m

This script will load the simulation scenarios, execute the BPbasedSLAM_dAoA_v2 algorithm, and get the estimated results.

---

## License
This project is licensed under the MIT License - see the LICENSE file for details.

