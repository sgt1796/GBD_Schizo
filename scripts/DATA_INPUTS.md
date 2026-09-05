# Data inputs

The study uses IHME Global Burden of Disease 2023 estimates for China and the
United States, 1990--2023.

## Canonical analysis files

- `GBD_data/GBD_2023_schizophrenia_fine_age_China_US.csv`: incidence,
  prevalence, and DALYs by year, country, sex, and five-year age group.
- `GBD_data/GBD_2023_population_China_US.csv`: matching population estimates
  used for count reconstruction and decomposition.

## Retained source archives

- `IHME-GBD_2023_DATA-774041bd-1.zip`: main fine-age schizophrenia export.
- `IHME-GBD_2023_DATA-22ef74c2-1.zip`: 70--74-year correction export.
- `IHME-GBD_2023_DATA-a9a792bb-1.zip`: population export.

All archives and canonical files are stored together in `GBD_data/`. The two
preparation scripts rebuild the canonical CSV files directly from the archives.

The source export contains reported zero values for incidence at ages 0--9 and
80 years or older, and for prevalence and DALYs at ages 0--9. These values are
retained as supplied. The primary decomposition uses the complete consecutive
five-year age range available for all required measures; the incidence-only
supported-age result is reported as a sensitivity analysis.
