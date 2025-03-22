# BondMatcher: H-Bond Stability Analysis in Molecular Systems

This repository contains the code and instructions for downloading the data used for the manuscript *"BondMatcher: H-Bond Stability Analysis in Molecular Systems"*.

## Summary

- [Installation](#installation)

- [Example](#running-an-example)

- [Database](#downloading-the-database)

## Installation

With a fresh install of Ubuntu 24.04 LTS

**1)** Download this repository  or run
   
```
 git clone https://github.com/thom-dani/BondMatcher.git
```

**2)** At the root of the repository, run

```
chmod +x install.sh
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
Next, enter the following command:
```
paraview example.pvsm
```
This example:
1. loads a database of extremum graphs
2. computes partial isomorphisms between them as well as the bond occurrence rate
3. displays the result in ParaView (see the two tabs `Layout #1` and `Layout #2`).
This reproduces the images from the Figure 8 of the manuscript *"BondMatcher: H-Bond Stability Analysis in Molecular Systems"*.

## Downloading the database
