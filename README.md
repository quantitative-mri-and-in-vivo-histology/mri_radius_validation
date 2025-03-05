# mri_radius_validation

This repository contains the code for the manuscript  
**"Towards MRI axon radius mapping in clinical settings: insights from MRI-scale histology and experimental validation".**

## Installation
The code requires Python 3.x and MATLAB. Additionally, please ensure the following software is installed:

### Python
```bash
pip install -r requirements.txt
```

Additionally, install **MATLAB Engine for Python**. Since the version must match your MATLAB installation, it must be installed manually. Follow the instructions at:  
[MATLAB Engine for Python](https://pypi.org/project/matlabengine)

### MATLAB
Ensure that the following toolboxes are installed and properly set up:
- **MISST Toolbox** (dMRI simulations)
- **NODDI Toolbox** (ex-vivo dMRI processing)
- **matlab_nifti** (ex-vivo dMRI processing)

### External Tools
These external tools must be installed and accessible in your system's `PATH`:
- **FSL** (dMRI processing) - [FSL website](https://fsl.fmrib.ox.ac.uk/fsl)
- **MRtrix3** (dMRI processing) - [MRtrix3 website](https://www.mrtrix.org)
- **ANTs** (dMRI processing) - [ANTs GitHub](https://github.com/ANTsX/ANTs)
- **c3d_affine_tool** (in-vivo dMRI processing) - [ITK-SNAP website](http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.C3D)

### Versions
The code has been tested with the following versions:
- FSL 6.0.6.4
- mrtix3 3.0.4
- ants 2.4.3
- MATLAB R2023b Update 8
- MISST Toolbox 0.93
- MODDI Toolbox 1.05
- python 3.10.16
- antspynet 0.2.9
- gradunwarp 1.2.3
- psutil 6.1.1
- matlabengine 23.2.3

## Hardware requirements
Most code can be run with standard hardware and 16-32 GB RAM. However, in-vivo processing requires a CUDA-capable GPU. Protocol optimization is computationally expensive and benefits from using a large number of cores typical for server CPUs. Runtime on a 72 core CPU was around 48 hours.

---
## Dataset
The dataset will be made available upon publication. The code assumes that the environment variable `MRV_DATA_PATH` points to the dataset root:

```bash
export MRV_DATA_PATH=/path/to/data
```
---

## Usage

### Processing of ex-vivo dMRI data
```bash
cd code
python -m dwi_processing.process_ex_vivo_dataset.py ${MRV_DATA_PATH}/mri_ex_vivo/processed
```

### Processing of in-vivo dMRI data
```bash
cd code
python -m dwi_processing.process_in_vivo_dataset.py ${MRV_DATA_PATH}/mri_in_vivo/processed
python -m dwi_processing.register_in_vivo_dataset.py ${MRV_DATA_PATH}/mri_in_vivo/processed
```

### Analyses presented in the manuscript
The code contains all scripts to reproduce the analyses and figures presented in the manuscript, located under `code/figures`. 

Some figures rely on simulations results shipped with the dataset under `${MRV_DATA_PATH}/simulations`. This data requires expensive computations (signal tables, protocol optimization) but can be regenerated running the following scripts:
```matlab
compute_signals_in_vivo_experimental.m   % Creates signal LUT for in-vivo dMRI simulations of experimental protocol
compute_signals_ex_vivo_experimental.m   % Creates signal LUT for ex-vivo dMRI simulations of experimental protocol
protocol_optimization_per_radius_simulate.m  % Creates signal LUTs for in-vivo protocol optimization
protocol_optimization_rois_simulate.m  % Mimics experimental radius estimation for in-vivo protocol candidates
```
---

## License
This project is licensed under the **MIT License**. However, some files or code snippets are subject to different licenses:
- **MPL 2.0**: `code/simulation/estimate_r_eff_non_linear_fit.m`
- **BSD-3-Clause**: `code/simulation/compute_lebedev_sphere.m`
- **Non-commercial research-only license** (New York University): `code/simulation/compute_even_order_spherical_harmonics.m`
