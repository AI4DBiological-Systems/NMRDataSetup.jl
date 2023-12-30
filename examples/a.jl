
#import CondaPkg # add if you're missing this from the Global Julia environment.
#import PythonCall # add if you're missing this from the Global Julia environment.

import PublicationDatasets as DS

import Random
Random.seed!(25)

using LinearAlgebra
using FFTW
using Printf

import PythonPlot
const PLT = PythonPlot

using Revise # comment this out if you don't have this installed.
import NMRDataSetup
const DSU = NMRDataSetup