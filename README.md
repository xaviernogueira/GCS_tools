# GCS Stage Analysis toolset
This repository contains a set of tools for processing geospatial data. While some tools may have a variety of fluvial geomorphological applications, an emphasis is placed on analysis of [geomorphic covariance structures](http://pasternack.ucdavis.edu/research/projects/geomorphic-covariance-structures/) (GCS). To launch the entire toolkit, run `master.py`. Otherwise, there is an individual GUI included for each tool.

#Prerequsites
-ArcGIS + Python 3.6 ; i.e. the `arcpy` python package and ArcGIS license.

-LASTools: functionality used for LiDAR Data Processing tool. LASTools can be downloaded [here](https://rapidlasso.com/lastools/).

## Getting Started -> two options

After installing any missing prerequisites, and downloading this repository you can either use:
 - 2D modeling and launch `master.py`, a description of the included tools are below (scripted by Kenneth Larrieu in Python 2.7 https://github.com/klarrieu)
   
   Note: The River GCS toolkit by Kenneth Larrieu is included as the GCS Analysis toolset builds off some of his functions. 
   If using the Python 2.7 toolset is desired, files labeled _oldpy are unapapted, but downloading from https://github.com/klarrieu/GCS_Scripts is recomended
   
 - Or follow the methodology outlined below to conduct a stage based GCS analysis (no 2D modeling required)
 
 ## Data requirements
 LiDAR data in LAZ or LAS format. Downloading the metadata includes shapefiles, which are necessary for defining the coorindate system or the LAS data
 
 Optional for precise vegetation classification based LiDAR processing:
    NAIP, 1m, four band imagery covering the LiDAR extent
    A generous (past floodplain) buffer of the reach centerline or a polygon cropped around the reach to avoid over processing ground points distal from the studied reach
 
 
