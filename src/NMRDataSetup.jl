__precompile__() # this module is safe to precompile
module NMRDataSetup

using LinearAlgebra
using FFTW
#using Statistics

import ForwardDiff
import Optim

# constant values.
function twopi(::Type{Float32})::Float32
    return 6.2831855f0 #convert(T, 2*π)
end

function twopi(::Type{Float64})::Float64
    return 6.283185307179586 #convert(T, 2*π)
end

function twopi(::Type{T})::T where T <: AbstractFloat
    return convert(T, 2*π)
end


include("DSP.jl")

include("bruker_IO.jl")

include("fit_singlet.jl")

include("utils.jl")

include("frontend.jl")

export 

    # General data containers.
    Data1D,
    OutputContainer1D,

    # data spectrum for fitting.
    SpectrumData1D,
    getSpectrumData1D,

    # fit singlet.
    FitSingletConfig,
    fitsinglet,

    # Bruker IO.
    
    Bruker1D1HSettings,
    loadBruker,
    setupBruker1Dspectrum,

    # utilities
    unpackcontainer,
    extractfreqrefHz,
    extractsolventfHz,
    getinversemap
end

