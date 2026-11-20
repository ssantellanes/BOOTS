---
title:'Codes and Supplemental Info: Santellanes et al. (2026) I hope Info'
date: 2026-11-20
permalink: README.md
---

# Santellanes and Melgar (2024) Codes
## Python Codes
This GitHub repo contains 4 python scripts. All except read_threadds_catalog.py are self contained and require some manipulation of the file paths.
In order to run read_threadds_catalog.py, do the following steps:

1. $ ./read_thredds_catalog.py -u https://www.ngdc.noaa.gov/thredds/catalog/dart_bpr/processed/catalog.xml -t nc > processed_dart_filelist_20200310.txt
2. $ wget --input-file processed_dart_filelist_20200310.txt [Note: delete the first and last lines of the file]

read_threadds_catalog.py was graciously provided by Aaron Sweeney from NOAA.

## Table S1

Table S1 contains all of the major data that was used in the manuscript. Depth is the latest depth from the DART station page. ppsd_serial.py extracts information for each DART station
for the depth at the day of the data sampling.

## Tsunami Database

The tsunami database from the National Geophysical Data Center / World Data Service: NCEI/WDS Global Historical Tsunami Database used for this study is attached here. 
See manuscript for which conditions were used to arrive at the tsunamis filtered for this study.

## Bayesian Denoiser

In order to run the Bayesian denoiser you will need to use the PPSDs which are in the NPZ folder. You will need to change the path to accomodate the paths, but other than that everything should work fine. The ipynb is pretty self-explanatory, so the code comments are kept to a minimum. The ERA5 files are needed to run anything, but they are nice to have around.

## MISC files
Files and uncritical files are included in this repo. They can be messed around with to include the ESNZ DART network or the SMART cable OBPGs. It's left open-ended here because my postdoc is ending, so maybe useful for someone to draft a paper (so long as they include me of course).
