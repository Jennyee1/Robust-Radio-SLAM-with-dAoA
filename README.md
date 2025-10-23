# Robust-Radio-SLAM-with-dAoA

This repository provides the official MATLAB implementation for the paper: \
**"ROBUST RADIO SLAM BY FUSING TEMPORAL-SPATIAL DIFFERENTIAL AOA INFORMATION COOPERATIVELY"** \
Submitted to: **ICASSP 2026**


---

## Requirements
MATLAB R2024b

## Usage

Current Status
This repository currently provides the core framework function described in the paper:
BPbasedSLAM_dAoA.m: This function implements the main Radio-SLAM framework algorithm, which fuses temporal and spatial dAoA information.

How to Run
A top-level main.m script or example file for reproducing the simulation results from the paper is not yet included. This is planned for a future update.

To use the code at this time, you will need to:
Prepare your own dataset (either simulated or experimental) according to the data structure required by the function.
Write your own main script to load the data and call the BPbasedSLAM_dAoA.m function to run the SLAM process.
We are working on providing a complete example script soon.

## License
This project is licensed under the MIT License.

