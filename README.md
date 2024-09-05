# SCHISM Input File Generation Scripts for UFS Coastal Application
## Project Structure
```
UFS-Coastal-Inputs/
│
├── data/
│   ├── hgrid.gr3
│   └── vgrid.in
│
├── outputs/
│   ├── albedo.gr3
│   ├── bctides.in
│   ├── diffmax.gr3
│   ├── diffmin.gr3
│   ├── manning.gr3
│   ├── watertype.gr3
│   └── windrot_geo2proj.gr3
│
├── scripts/
│   ├── coastal_ike_shinnecock_atm2sch/
│   │   ├── gen_bctides.py
│   │   ├── gen_bnd.py
│   │   ├── gen_gr3_input.py
│   │   └── README
│   │
│   └── Test_Duck_NC/
│       └── README
│
├── LICENSE
└── README.md
```

This repository contains Python scripts for generating input files used in SCHISM for the UFS Coastal regression tests.

## Available Scripts

### 1. Coastal Ike Shinnecock

- **Location:** `/scripts/coastal_ike_shinnecock_atm2sch`

### 2. Duck, North Carolina Test
- **Location:** `/scripts/Test_Duck_NC`



