#!/usr/bin/env python

import sys

if len(sys.argv) >= 1:
    inputPath = sys.argv[1]
else:
    print("Usage: ./separatrixExtraction.py <input file>")
    sys.exit()

# state file generated using paraview version 5.13.0
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------
#start_datavti = XMLImageDataReader(FileName=['start_data.vti'])
start_datavti = XMLPartitionedImageDataReader(FileName=[inputPath])

## create a new 'Pass Arrays'
#passArrays1 = PassArrays(registrationName='PassArrays1', Input=start_datavti)
#passArrays1.PointDataArrays = ['rho']

## create a new 'Calculator'
#calculator1 = Calculator(registrationName='Calculator1', Input=passArrays1)
#calculator1.ResultArrayName = 'opposite'
#calculator1.Function = '-rho'
## WARNING: use Double here! (otherwise diagonals)
#calculator1.ResultArrayType = 'Double'

## create a new 'Pass Arrays'
#passArrays2 = PassArrays(registrationName='PassArrays2', Input=calculator1)
#passArrays2.PointDataArrays = ['opposite']

## create a new 'Threshold'
#threshold1 = Threshold(registrationName='Threshold1', Input=passArrays2)
#threshold1.Scalars = ['POINTS', 'opposite']
#threshold1.LowerThreshold = -157.53614805600088
#threshold1.UpperThreshold = -0.01

## create a new 'Resample To Image'
#resampleToImage1 = ResampleToImage(registrationName='ResampleToImage1', Input=threshold1)
#resampleToImage1.SamplingDimensions = [int(resolution), int(resolution), int(resolution)]
#resampleToImage1.SamplingBounds = [16.08506202697754, 23.57298469543457, -42.79999542236328, -35.95000457763672, -20.957918167114258, -14.282181739807129]

## create a new 'Resample With Dataset'
#resampleWithDataset1 = ResampleWithDataset(registrationName='ResampleWithDataset1', SourceDataArrays=calculator1,
#    DestinationMesh=resampleToImage1)
#resampleWithDataset1.PassCellArrays = 1
#resampleWithDataset1.PassPointArrays = 1
#resampleWithDataset1.PassPartialArrays = 1
#resampleWithDataset1.MarkBlankPointsAndCells = 1
#resampleWithDataset1.SnapToCellWithClosestPoint = 1
#resampleWithDataset1.CellLocator = 'Static Cell Locator'

tTKTopologicalSimplificationByPersistence1 = TTKTopologicalSimplificationByPersistence(Input=start_datavti)
tTKTopologicalSimplificationByPersistence1.InputArray = ['POINTS', 'opposite']
tTKTopologicalSimplificationByPersistence1.PersistenceThreshold = 1e-03

# create a new 'TTK MorseSmaleComplex'
tTKMorseSmaleComplex1 = TTKMorseSmaleComplex(registrationName='TTKMorseSmaleComplex1', Input=tTKTopologicalSimplificationByPersistence1)
tTKMorseSmaleComplex1.ScalarField = ['POINTS', 'opposite']
tTKMorseSmaleComplex1.OffsetField = ['POINTS', 'opposite']
tTKMorseSmaleComplex1.DiscreteGradientBackend = 'Stochastic Discrete Gradient'
tTKMorseSmaleComplex1.AscendingSegmentation = 0
tTKMorseSmaleComplex1.DescendingSegmentation = 0
tTKMorseSmaleComplex1.MorseSmaleComplexSegmentation = 0
tTKMorseSmaleComplex1.ReturnSaddleConnectors = 1
tTKMorseSmaleComplex1.SaddleConnectorsPersistenceThreshold = 1e-04

# create a new 'Threshold'
threshold3 = Threshold(Input=OutputPort(tTKMorseSmaleComplex1, 1))
threshold3.Scalars = ['CELLS', 'SeparatrixType']
threshold3.LowerThreshold = 0
threshold3.UpperThreshold = 0

SaveData("start_data.vtu", threshold3)
