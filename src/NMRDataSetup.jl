module NMRDataSetup

import PyCall, Optim, NLopt

include("utils.jl")
include("DSP.jl")
include("load.jl")
include("assemble_mixture.jl")

export loadspectrum, getwraparoundDFTfreqs, evalcomplexLorentzian, gettimerange

end
