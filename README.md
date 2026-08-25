<div align="center">
 <h1>metaDRMetRap</h1>
</div>

<!-- badges: start -->
[![Paper](https://img.shields.io/badge/paper-10.1111%2Facel.70131-blue)](https://onlinelibrary.wiley.com/doi/10.1111/acel.70131)
[![Data & Code](https://img.shields.io/badge/Zenodo-15673918-blue)](https://zenodo.org/records/15673918)
<!-- badges: end -->

Data metadata for **"Rapamycin, not metformin, mirrors dietary restriction-driven
lifespan extension in vertebrates: a meta-analysis"**.

- Paper: <https://onlinelibrary.wiley.com/doi/10.1111/acel.70131>
- Data and code: <https://zenodo.org/records/15673918>

## Overview

Dietary restriction reliably extends lifespan across a wide range of species, and
two drugs — rapamycin and metformin — are widely promoted as pharmacological
mimetics of that effect. Whether either actually reproduces the DR response in
vertebrates had not been assessed systematically.

This meta-analysis compares lifespan extension under dietary restriction,
rapamycin, and metformin across vertebrate studies, drawing on 911 effect sizes.
It finds that rapamycin, but not metformin, mirrors the DR-driven extension of
lifespan. This repository holds the **data dictionary** for that synthesis; the
full data and code are archived on Zenodo.

## Authors

Edward R. Ivimey-Cook\*@, Zahida Sultanova\*@, and Alexei A. Maklakov

\*these authors contributed equally
@corresponding authors: <e.ivimeycook@gmail.com> and <zahida.sultanova@uea.ac.uk>

## Files

| File | Contents |
| :--- | :------- |
| `analysis_data.csv` | The meta-analytic dataset (911 effect sizes) |
| `Analysis_Script.R` | Master analysis script; all other scripts are sourced within this one |
| `session.txt` | R version and packages loaded, with versions |

Create a `Code` folder and run the script `Analysis_Script.R`.

## Data dictionary

`analysis_data.csv` — all column names are cleaned when loading.

**Study identification**

| Column | Description |
| :----- | :---------- |
| `id` | 1:911 unique individual ID |
| `author list` | List of authors |
| `title` | Title of paper |
| `Year` | Year paper was published (extracted from Scopus/WoS) |

**Control group**

| Column | Description |
| :----- | :---------- |
| `m_control` | Control lifespan measure |
| `n_control` | Sample size for control |
| `sd_control` | Standard deviation for control (left blank if not available) |
| `se_control` | Standard error for control (left blank if not available) |

**Treatment group**

| Column | Description |
| :----- | :---------- |
| `m_treatment` | Treatment lifespan measure |
| `n_treatment` | Sample size for treatment |
| `sd_treatment` | Standard deviation for treatment (left blank if not available) |
| `se_treatment` | Standard error for treatment (left blank if not available) |

**Measurement details**

| Column | Description |
| :----- | :---------- |
| `m_Measure` | Median or mean lifespan |
| `m_raw` | Whether raw data was given to calculate measures |
| `m_time` | Unit of lifespan measure |
| `m_Location` | Location of lifespan measure in the paper |

**Study characteristics**

| Column | Description |
| :----- | :---------- |
| `Species` | Species sampled |
| `Treatment_Type` | Type of treatment used (e.g. fasting, percent reduction, rapamycin) |
| `Treatment_Overall` | Broad treatment groups (Rapamycin, Metformin, DR) |
| `Treatment_Name` | Treatment name in paper |
| `Sex` | Sex studied |
| `Strain` | The name of any strain used, if appropriate |
| `Other Variables` | Whether any other variables were tested in the experiment |
| `Notes` | Additional notes |
| `mice900_keep` | Whether or not the mice paper passed the 900 day rule |

## R environment

R version and packages loaded (with versions): see the `session.txt` file.

## Citation

> Ivimey-Cook, E. R., Sultanova, Z., & Maklakov, A. A. (2025). Rapamycin, not
> metformin, mirrors dietary restriction-driven lifespan extension in
> vertebrates: a meta-analysis. *Aging Cell*.
> <https://doi.org/10.1111/acel.70131>

## Contact

Edward R. Ivimey-Cook — <e.ivimeycook@gmail.com> —
[ORCID 0000-0003-4910-0443](https://orcid.org/0000-0003-4910-0443)
