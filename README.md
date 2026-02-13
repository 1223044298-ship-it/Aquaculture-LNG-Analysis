# LNG Cold Energy Unlocks Climate-Resilient Mariculture for Low-Carbon Nutrition

This repository contains the code and processed datasets used in the manuscript:

**"LNG cold energy unlocks climate-resilient mariculture for low-carbon nutrition"**

Submitted to *Nature Food*.

All analyses are fully reproducible using the scripts and publicly available datasets described below. No proprietary software was used.

---

## Overview

This repository provides:

- National-scale scenario modeling scripts
- Thermal suitability assessment workflows (~10 m upper-ocean temperature)
- Essential amino acid (EAA) supply calculations
- Self-sufficiency ratio (SSR) computations
- Child stunting reduction estimations
- Carbon and economic performance calculations
- Figure-generation scripts (Figures 1–3)

The study links LNG regasification cold discharge to potential nearshore aquaculture expansion, improved essential amino acid supply, enhanced national self-sufficiency, and reduced carbon intensity of nutrition pathways.

---

## Repository Structure

```
Aquaculture-LNG-Analysis/
│
├── data/
│   ├── input/
│   │   ├── SSRcapped_output.xlsx
│   │   ├── child.xlsx
│   │   ├── country_eaa_weighted_recommendations_and_coverage_with_total.xlsx
│   │   ├── site_species_cover10_needcool.xlsx
│   │   ├── species_temp_ranges.xlsx
│   │
│   ├── shapefiles/
│   │   └── ne_50m_admin_0_countries/
│
├── scripts/
│   ├── fig1_global_map.py
│   ├── fig2_complete_maps.py
│   ├── fig2_radial_panels.R
│   ├── fig3_left_iconpanel.py
│   ├── fig3_right_panel.py
│
└── README.md
```

---

## Data Sources

All datasets used in this study are publicly available:

- CMEMS global ocean reanalysis (monthly upper-ocean temperature, ~10 m depth)
- FAO FishStatJ production and trade statistics
- WHO/FAO amino acid requirement guidelines
- United Nations population datasets
- UNICEF/World Bank child stunting statistics
- Global LNG terminal infrastructure databases

Raw NetCDF temperature files are not included in this repository due to size constraints but can be downloaded directly from the Copernicus Marine Environment Monitoring Service (CMEMS) portal.

Processed datasets required to reproduce figures are included in the `data/` directory.

---

## System Requirements

### Software

- Python 3.9 or later
- R 4.2 or later

### Python Packages

- numpy
- pandas
- matplotlib
- geopandas
- shapely
- xarray
- cartopy

### R Packages

- tidyverse
- readxl
- gghalves
- ggplot2

All scripts were tested on:

- Windows 10/11

No non-standard hardware is required.

---

## Installation

1. Clone the repository:

```bash
git clone https://github.com/1223044298-ship-it/Aquaculture-LNG-Analysis.git
```
2. Navigate into the repository directory:

```bash
cd Aquaculture-LNG-Analysis
```

3. (Optional but recommended) Create a virtual environment:

```bash
python -m venv venv
venv\Scripts\activate     # Windows
```

4. Install required Python packages:

```bash
pip install numpy pandas matplotlib geopandas shapely xarray cartopy
```

5. Install required R packages (inside R console):

```R
install.packages(c("tidyverse","readxl","gghalves","ggplot2"))
```

Typical installation time on a standard desktop computer:  
**< 10 minutes**

---

## Running the Code

### Figure 1 – Global Thermal Suitability Map

```bash
python scripts/fig1_global_map.py
```

### Figure 2 – Global EAA Coverage Maps

```bash
python scripts/fig2_complete_maps.py
```

### Figure 2 Radial Panels

```bash
Rscript scripts/fig2_radial_panels.R
```

### Figure 3 Panels

```bash
python scripts/fig3_left_iconpanel.py
python scripts/fig3_right_panel.py
```

Expected runtime per script on a standard desktop computer:  
**1–5 minutes**

---

## Expected Output

Running the scripts will generate:

- Global thermal suitability maps
- National EAA coverage maps
- Radial EAA composition panels
- Child stunting reduction visualizations

Outputs are saved locally in the working directory unless otherwise specified in each script.

---

## Reproducibility

All model assumptions, equations, and scenario parameters are documented in the manuscript Methods section.  

All analyses are deterministic and rely exclusively on publicly available datasets.  

No proprietary software or closed-source tools were used.

The repository includes all processed datasets required to reproduce manuscript figures.

---

## Citation

If using this repository, please cite:

> [Authors]. LNG cold energy unlocks climate-resilient mariculture for low-carbon nutrition. *Nature Food* (under review).

---

## License

This repository is made available for academic and research purposes.  
Please contact the corresponding author for additional permissions if required.
