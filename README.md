# Dissertation Code and Results

This repository contains the R code, input data, processed block maxima datasets, model outputs, diagnostic plots, and return level results used for my B.Sc. (Hons.) dissertation:

**Non-Stationary Extreme Value Analysis of Climate Extremes in the Netherlands**

## Overview

The dissertation applies Extreme Value Theory (EVT) to Dutch climate extremes using the Block Maxima approach. Stationary and non-stationary Generalised Extreme Value (GEV) models are fitted and compared for climate variables recorded at two Royal Netherlands Meteorological Institute (KNMI) stations:

- **De Bilt**: inland station
- **De Kooy**: coastal station

The variables analysed are:

- Annual maximum daily temperature
- Annual minimum daily temperature
- Annual maximum daily precipitation
- Annual maximum hourly wind speed

For each dataset, a stationary GEV model and seven non-stationary GEV specifications are fitted. Model comparison is carried out using negative log-likelihood, AIC, BIC, and likelihood ratio tests. Diagnostic plots and return level estimates are also produced.

## Repository Structure

```text
Dissertation-Code-and-All-Results/
│
├── Codes/
│   └── R scripts used for data processing, model fitting, diagnostics, and return level analysis.
│
├── raw knmi obs/
│   └── Raw KNMI observational data used in the analysis.
│
├── bm_output_debilt_tmax_FIXED/
│   └── Output files for De Bilt annual maximum daily temperature.
│
├── bm_output_debilt_minima_FIXED/
│   └── Output files for De Bilt annual minimum daily temperature.
│
├── bm_output_debilt_precip_FIXED/
│   └── Output files for De Bilt annual maximum daily precipitation.
│
├── bm_output_debilt_wind_FIXED/
│   └── Output files for De Bilt annual maximum hourly wind speed.
│
├── bm_output_dekooy_tmax_FIXED/
│   └── Output files for De Kooy annual maximum daily temperature.
│
├── bm_output_dekooy_tmin_FIXED/
│   └── Output files for De Kooy annual minimum daily temperature.
│
├── bm_output_dekooy_precip_FIXED/
│   └── Output files for De Kooy annual maximum daily precipitation.
│
├── bm_output_dekooy_wind_fhx_FIXED/
│   └── Output files for De Kooy annual maximum hourly wind speed.
│
└── README.md
