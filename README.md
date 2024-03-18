---
title:'The 2020 Sand Point Earthquake: (Disaster) Chart Party'
date: 2024-03-18
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
