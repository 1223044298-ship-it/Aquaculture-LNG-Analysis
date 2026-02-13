# LNG Cold Energy Unlocks Climate-Resilient Mariculture for Low-Carbon Nutrition

This repository contains the code and processed datasets used in the manuscript:

"LNG cold energy unlocks climate-resilient mariculture for low-carbon nutrition"

Submitted to *Nature Food*.

All analyses are fully reproducible using the scripts and publicly available datasets described below. No proprietary software was used.

---

## 1. Overview

This repository provides:

- National-scale scenario modeling scripts
- Thermal suitability assessment workflows (~10 m upper-ocean temperature)
- Essential amino acid (EAA) supply calculations
- Self-sufficiency ratio (SSR) computations
- Child stunting reduction estimations
- Carbon and economic performance calculations
- Figure-generation scripts (Figures 1–3)

The study links LNG regasification cold discharge to potential nearshore aquaculture expansion, improved essential amino acid supply, enhanced national self-sufficiency, and reduced carbon intensity of nutrition pathways.

All scripts are deterministic and generate identical outputs when run with the included datasets.

---

## 2. Repository Structure

Aquaculture-LNG-Analysis/

- data/
  - input/
    - SSRcapped_output.xlsx
    - child.xlsx
    - country_eaa_weighted_recommendations_and_coverage_with_total.xlsx
    - site_species_cover10_needcool.xlsx
    - species_temp_ranges.xlsx
  - shapefiles/
    - ne_50m_admin_0_countries/

- scripts/
  - fig1_global_map.py
  - fig2_complete_maps.py
  - fig2_radial_panels.R
  - fig3_left_iconpanel.py
  - fig3_right_panel.py

- README.md

---

## 3. Data Sources

All datasets used in this study are publicly available:

- Copernicus Marine Environment Monitoring Service (CMEMS) global ocean reanalysis (monthly upper-ocean temperature, ~10 m depth, variable: thetao)
- FAO FishStatJ production and trade statistics
- WHO/FAO amino acid requirement guidelines
- United Nations population datasets
- UNICEF/World Bank child stunting statistics
- Global LNG terminal infrastructure databases

### Large NetCDF Files

Raw CMEMS NetCDF temperature files are not included in this repository due to size constraints.  

They can be downloaded directly from:

https://marine.copernicus.eu/

The scripts are fully compatible with standard CMEMS monthly reanalysis products (variable: thetao, ~10 m depth).

Processed datasets required to reproduce all manuscript figures are included in the `data/` directory.

---

## 4. Demo Dataset

The repository includes small processed input datasets in:

data/input/

These Excel files are sufficient to:

- Run all figure-generation scripts
- Reproduce national EAA coverage calculations
- Reproduce SSR outputs
- Generate child stunting impact figures

No large raw data files are required to reproduce the published figures.

---

## 5. System Requirements

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

Tested on:

- Windows 10/11

No non-standard hardware is required.

---

## 6. Installation

1. Clone the repository:

git clone https://github.com/1223044298-ship-it/Aquaculture-LNG-Analysis.git

2. Navigate into the repository:

cd Aquaculture-LNG-Analysis

3. (Optional) Create virtual environment:

python -m venv venv
venv\Scripts\activate   # Windows

4. Install Python packages:

pip install numpy pandas matplotlib geopandas shapely xarray cartopy

5. Install R packages (inside R console):

install.packages(c("tidyverse","readxl","gghalves","ggplot2"))

Typical installation time on a standard desktop computer:
< 10 minutes

---

## 7. Running the Code

### Figure 1 – Global Thermal Suitability Map

python scripts/fig1_global_map.py

### Figure 2 – Global EAA Coverage Maps

python scripts/fig2_complete_maps.py

### Figure 2 – Radial EAA Panels

Rscript scripts/fig2_radial_panels.R

### Figure 3 – Child Stunting Panels

python scripts/fig3_left_iconpanel.py
python scripts/fig3_right_panel.py

Expected runtime per script on a standard desktop computer:
1–5 minutes

---

## 8. Expected Output

Running the scripts will generate:

- Global thermal suitability maps
- National EAA coverage maps
- Radial EAA composition panels
- Child stunting reduction visualizations
- Scenario comparison graphics

Outputs are saved locally in the working directory unless otherwise specified in the script.

---

## 9. Reproducibility

All model assumptions, equations, and scenario parameters are fully documented in the manuscript Methods section.

All analyses are deterministic.

All processed datasets required to reproduce manuscript figures are included.

Raw source datasets are publicly available from the providers listed above.

No proprietary software or closed-source tools were used.

---

## 10. License

This repository is provided for academic and research purposes.

Please cite the associated manuscript if using this code in derivative work.
