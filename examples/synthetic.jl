### load a few 700 MHz BMRB experiments, each with a known concentration.
# combine them into a single FID.

import NMRData


#include("../src/NMRDataSetup.jl")
import NMRDataSetup

using FFTW
import PyPlot
import BSON

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

# user inputs.
BMRB_setup_label = "BMRB-700" # for use with the NMRData package. See that package's readme for details.

save_dir = "/home/roy/MEGAsync/outputs/NMR/combined/"
project_name = "cal01"

target_names = ["D-(+)-Glucose"; # compounds you wish to add to the mixture.
"L-Alanine"]


weights = ones(Float64, length(target_names)) # mixing weights for each experiment.

solvent_ppm_guess = 4.7 # guess where the solvent ppm is for each BMRB experiment.
solvent_window_ppm = 0.1 # allowed solvent shift, in ppm.

### end inputs.




## get BMRB experiment paths from the NMRData package. You can specify your own array of experiments to mix.
paths = Vector{String}(undef, 0)
for i = 1:length(target_names)
    push!( paths, NMRData.getBMRBexperimentpath(BMRB_setup_label, target_names[i])[1] )
end
@assert all(isdir.(paths)) # check for bad paths.


## mix.
s_t, ν_0ppm, fs, SW,
s_t_set, ν_0ppm_set, fs_set, SW_set = NMRDataSetup.assemblemixture(paths; weights = weights)


## store.
save_folder_path = joinpath(save_dir, project_name)
isdir(save_folder_path) || mkdir(save_folder_path); # make save folder if it doesn't exist.

save_path = joinpath(save_folder_path, "$(project_name).bson")
BSON.bson(save_path,
s_t = s_t,
molecule_names = target_names,
fs = fs,
SW = SW,
ν_0ppm = ν_0ppm,
s_t_set = s_t_set,
ν_0ppm_set = ν_0ppm_set,
fs_set = fs_set,
SW_set = SW_set)


#### plot.

function getDFTfreqrange(N::Int, fs::T)::LinRange{T} where T
    a = zero(T)
    b = fs-fs/N

    return LinRange(a, b, N)
end

# Assumes U0 is sorted in ascending order.
function wrapfreqrange(U0, ν_begin::T, fs::T) where T <: Real

    N = length(U0)
    ind = findfirst(xx->(xx>ν_begin), U0)

    # TODO handle these exceptions with more grace.
    @assert typeof(ind) == Int
    @assert ind <= N

    out = Vector{T}(undef, N)
    #out[1:ind] = U0[1:ind] .+ fs
    #out[ind+1:end] = U0[ind+1:end]

    M = N-ind
    out[1:M] = U0[ind+1:end]
    out[M+1:end] = U0[1:ind] .+ fs

    inds = collect(1:N)
    inds[1:M] = collect(ind+1:N)
    inds[M+1:end] = collect(1:ind)

    return out, inds
end

function getwraparoundDFTfreqs(N::Int, fs::T, ν_begin::T) where T

    U0 = getDFTfreqrange(N, fs)
    out, inds = wrapfreqrange(U0, ν_begin, fs)

    return U0, out, inds
end

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

offset_Hz = ν_0ppm - (ppm2hzfunc(0.3)-ppm2hzfunc(0.0))

N = length(s_t)
DFT_s = fft(s_t)
U_DFT, U, U_inds = getwraparoundDFTfreqs(N, fs, offset_Hz)
y = DFT_s[U_inds]

##############################

S_U = DFT_s[U_inds]
P = hz2ppmfunc.(U)

S_set_U = collect( (fft(s_t_set[n]) ./fs)[U_inds] for n = 1:length(s_t_set))
#P_cost = hz2ppmfunc.(U_cost)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(S_U), label = "combined spectrum")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("simulated spectra")


