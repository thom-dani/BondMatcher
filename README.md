# BondMatcher: H-Bond Stability Analysis in Molecular Systems

This repository contains the code and instructions for downloading the data used for the manuscript *"BondMatcher: H-Bond Stability Analysis in Molecular Systems"*.

These instructions have been tested with a fresh install of Ubuntu 24.04 LTS. Adjustments may be required for other operating systems.

## Summary

- [Installation](#installation)

- [Example](#running-an-example)

- [Database](#downloading-the-database)

## Installation

**1)** Download this repository  or run
   
```
 git clone https://github.com/thom-dani/BondMatcher.git
```

**2)** At the root of the repository, run

```
./install.sh
```

**3) Optional** : To permanently add the paraview and ttk environment variables to your .bashrc, run
```
echo 'PV_PREFIX=$(pwd)/ttk-paraview/install' >> ~/.bashrc
echo 'export PATH=$PATH:$PV_PREFIX/bin' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=$PV_PREFIX/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
echo 'export PYTHONPATH=$PYTHONPATH:$PV_PREFIX/lib/python3.12/site-packages' >> ~/.bashrc
```
and 
```
echo 'TTK_PREFIX=$(pwd)/ttk/install' >> ~/.bashrc
echo 'export PV_PLUGIN_PATH=$TTK_PREFIX/bin/plugins/TopologyToolKit' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=$TTK_PREFIX/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
echo 'export PYTHONPATH=$PYTHONPATH:$TTK_PREFIX/lib/python3.12/site-packages' >> ~/.bashrc
```
at the root of this repository and 
```
source ~/.bashrc
```
to apply the changes.


**4)** At this point, you can launch paraview from the command line with 

```
paraview
```
or use the following scripts to automatically run the examples.

## Running an example
To run the provided example, first, go to the `example` directory (this is important as all paths are stored in the example file relatively to that directory).
Next, from there, enter the following command:
```
paraview example.pvsm
```
This example:
1. loads a database of extremum graphs
2. computes partial isomorphisms between them as well as the bond occurrence rate
3. displays the result in ParaView (see the two tabs `Layout #1` and `Layout #2`).


This reproduces the images from the Figure 8 of the manuscript *"BondMatcher: H-Bond Stability Analysis in Molecular Systems"*. To render the extremum graphs with straight lines (as done in the manuscript), increase the number of iterations of the pipeline entry `TTKGeometrySmoother1` to a sufficiently high value (e.g., > 1,000).

## Downloading the database
We contribute a database of 4544 electron densities (computed on a `256**3` grid,
around 800 GB of decompressed data, see the companion manuscript for futher details). It is organized as follows:

- Two pathways of proton tunneling for the *Prism* water hexamer:
  - With anti-geared motion (256 steps)
    - [Download file (32.7 GB)](https://zenodo.org/records/14909099/files/W6_bifdrop_sp_along_tunneling_pathways.tar.gz)
    - To decompress the file and prepare the corresponding cinema database, move the downloaded file to the `scripts/pathways` directory. Next, go to the directory `scripts/pathways`. From there, enter the following command:
    ```
    ./prepareCinemaDataBase.sh W6_bifdrop_sp_along_tunneling_pathways.tar.gz
    ```
    This will create a cinema data base of electron densities (in *.vti format) which can be interactively explored within ParaView using the TTK Cinema filters (typically, `TTKCinemaReader` to read the database, followed by `TTKCinemaQuery` to select entries, followed by `TTKCinemaProductReader` to load the selected densities, see the above example).
    
    Optionally, to create the cinema database of corresponding extremum graphs, enter the following command (this will take a **LONG** time, if needed decrease the resolution parameter from `512` to lower values):
    ```
    ./separatrixExtraction.sh  W6_bifdrop_sp_along_tunneling_pathways.cdb 512
    ```
  - With geared motion (256 steps)
    - [Download file (32.7 GB)](https://zenodo.org/records/14909089/files/W6_drop1_sp_along_tunneling_pathways.tar.gz)
    - To decompress the file and prepare the corresponding cinema database, move the downloaded file to the `scripts/pathways` directory. Next, go to the directory `scripts/pathways`. From there, enter the following command:
    ```
    ./prepareCinemaDataBase.sh W6_drop1_sp_along_tunneling_pathways.tar.gz
    ```
    This will create a cinema data base of electron densities (in *.vti format) which can be interactively explored within ParaView using the TTK Cinema filters (typically, `TTKCinemaReader` to read the database, followed by `TTKCinemaQuery` to select entries, followed by `TTKCinemaProductReader` to load the selected densities, see the above example).
    
    Optionally, to create the cinema database of corresponding extremum graphs, enter the following command (this will take a **LONG** time, if needed decrease the resolution parameter from `512` to lower values):
    ```
    ./separatrixExtraction.sh  W6_drop1_sp_along_tunneling_pathways.cdb 512
    ```

- 48 vibrational perturbations for the *Ring*, *Book*, *Cage* and *Prism* water hexamers:
  - *Ring* hexamer:
    - [Download file 1 (42.9 GB)](https://zenodo.org/records/14911658/files/W6_ring_refined_DFT_ScaledRange4.0_Ngeoms21.tar.gz.partaa)
    - [Download file 2 (42.9 GB)](https://zenodo.org/records/14911654/files/W6_ring_refined_DFT_ScaledRange4.0_Ngeoms21.tar.gz.partab)
    - [Download file 3 (42.9 GB)](https://zenodo.org/records/14911650/files/W6_ring_refined_DFT_ScaledRange4.0_Ngeoms21.tar.gz.partac)
    - [Download file 4 (198.6 MB)](https://zenodo.org/records/14911646/files/W6_ring_refined_DFT_ScaledRange4.0_Ngeoms21.tar.gz.partad)
