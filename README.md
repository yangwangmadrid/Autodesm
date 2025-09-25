# Autodesm
Autodesmotic reactions for strain energy evaluation

**Source codes:** Autodesmotic reactions for strain energy evaluation  
**Last modified:** September 25, 2025  
**License:** For academic and non-commercial use only  
**Author:** Yang Wang (ORCID: 0000-0003-2540-2199; Email: yangwang@yzu.edu.cn)

---

## Overview

This repository provides the datasets, MATLAB source codes, and demonstration scripts used in the study:  

**"Autodesmotic reactions for general strain energy evaluation in polycyclic aromatic nanocarbons."**

The package includes:

- Datasets used for training the models
- MATLAB scripts for training models and applying them to new systems  
- Trained models for energy and bond length prediction  
- Benchmark test cases for carbon nanobelts (CNBs), carbon nanotubes (CNTs),
  helicenes, bowl-shaped PAHs, and fullerenes  

All molecular geometries of polycyclic aromatic hydrocarbons (PAHs) and test
systems were optimized at the B3LYP/6-31G(d) level, and all total energies were
obtained from B3LYP/6-311G(d)//B3LYP/6-31G(d) single-point calculations.  

All scripts have been tested and verified with **MATLAB R2023**.


---

## Repository Structure

```text
.
├── Datasets/                     # Molecular datasets
│   ├── Benzenoids_C48H24/        # 1,516 isomers of C48H24
│   │   ├── DATA_complete.dat     # Data for energy model training
│   │   └── AllData_BL_C48H24.mat # Data for bond length model training
│   └── Benzenoids_C6H6_C96H24/   # 2,275 PAHs from C6H6 to C96H24
│       ├── DATA_small.dat        # Data for energy model training
│       ├── DATA_medium.dat       # Data for energy model training
│       ├── DATA_larger.dat       # Data for energy model training
│       ├── DATA_big.dat          # Data for energy model training
│       ├── DATA_big88.dat        # Data for energy model training
│       └── AllData_BL.mat        # Data for bond length model training
│
├── Models_Train/                 # Model training scripts
│   ├── Energy_Models/
│   │   ├── Isomeric_C48H24/
│   │   │   └── model_EDFT_EHMO_CCBond_numHH.m
│   │   └── Heterogeneous/
│   │       └── model_EDFT_EHMO_CCBond_numHH.m
│   └── BondLength_Models/
│       ├── Isomeric_C48H24/
│       │   └── model_BL_fit.m
│       └── Heterogeneous/
│           └── model_BL_fit.m
│
└── Test_Cases/                   # Demonstration test cases
    ├── Fig3_Benchmark_CNBs_C48H24/
    │   ├── CNB_6_6_1.m
    │   ├── ... ... ...
    │   └── CNB_8_4_21.m
    ├── Fig4_Benchmark_Armchair_CNBs/
    │   └── armchair_CNBs.m
    ├── Fig5_CNBs_Helicenes_Bowls/
    │   ├── a_CNB_9_0_3_1.m
    │   ├── b_7circulene.m
    │   ├── ... ... ...
    │   └── d_carboncone.m
    ├── Fig6_CNTs/
    │   ├── a_cnt_5_5_2.m
    │   ├── a_cnt_5_5_21.m
    │   ├── ... ... ...
    │   └── b_cnt_12_6_1.m
    ├── Fig7_Fullerenes/
    │   ├── C60_2.m
    │   ├── ... ... ...
    │   └── C60_1812.m
    └── Vogtle_Belt/
        ├── vogtle_belt.m
        └── mobius_belt.m
```

---

## Datasets

### Isomeric Space of C48H24

- DATA_complete.dat
  ASCII file with information for 1,516 planar benzenoid PAHs of C48H24:
    - Encoded name (BEC code)
    - Planarity
    - DFT total energy (in Hartree)
    - Simple HMO energy (in |beta|)
    - Atom counts (C, H)
    - H···H count and distances (in Angstrom)
    - C–C bond count and bond lengths (in Angstrom)
    - C–H bond count and bond lengths (in Angstrom)
  
  Each line corresponds to each PAH molecule.

- AllData_BL_C48H24.mat
  MATLAB binary file containing bond-level variables for 90,960 C–C bonds:
    - `bndLen`: C–C bond lengths (in Angstrom)
    - `bndOrd`: Coulson bond orders (distance-dependent HMO)
    - `bndTyp`: Simple C–C bond types (CH–CH, CH–C, C–C)
    - `bndTypExt`: Atom-based bond types
    - `bndTypRB`: Ring-based bond types


### Heterogeneous Space (C6H6 to C96H24)

- Data files (in the same format as Benzenoids_C48H24/DATA_complete.dat):
    - `DATA_small.dat` (C6H6 – C34H20)
    - `DATA_medium.dat` (C36H18 – C56H24)
    - `DATA_larger.dat` (C56H24 – C82H24)
    - `DATA_big.dat` (C84H24 – C96H24)
    - `DATA_big88.dat` (C84H24, C88H26)
  
  A total of 2,275 planar benzenoid PAHs

- AllData_BL.mat
  
  MATLAB binary file containing bond-level variables for 205,350 C–C bonds
across the heterogeneous dataset.
  Format is identical to `AllData_BL_C48H24.mat`.


---

## Model Training

### Energy Models

  - Isomeric space of C48H24:
    
    Run `model_EDFT_EHMO_CCBond_numHH.m` in `Energy_Models/Isomeric_C48H24/` to train the energy prediction model.

  - Heterogeneous space (C6H6–C96H24):
    
    Run `model_EDFT_EHMO_CCBond_numHH.m` in `Energy_Models/Heterogeneous/`.

### Bond Length Models

  - Isomeric space of C48H24:
    
    Run `model_BL_fit.m` in `BondLength_Models/Isomeric_C48H24/` to train the bond length prediction model.

  - Heterogeneous space (C6H6–C96H24):
    
    Run `model_BL_fit.m` in `BondLength_Models/Heterogeneous/`.


---

## Test Cases

Benchmark calculations can be reproduced using the provided MATLAB scripts:
  - C48H24 CNBs (Fig. 3): `Test_Cases/Fig3_Benchmark_CNBs_C48H24/`
  - Armchair CNBs of varying size (Fig. 4): run `armchair_CNBs.m`
  - Unconventional CNBs, helicenes, and bowls (Fig. 5): `Test_Cases/Fig5_CNBs_Helicenes_Bowls/`
  - CNTs (Fig. 6): `Test_Cases/Fig6_CNTs/`
    
    Note: Only two representative lengths per CNT type are provided.
  - C60 fullerenes (Fig. 7): `Test_Cases/Fig7_Fullerenes/`
    
    Note: Selected representative isomers covering different NAPPs (0–16).
  - Vögtle belt and its Möbius form: `Test_Cases/Vogtle_Belt/`

