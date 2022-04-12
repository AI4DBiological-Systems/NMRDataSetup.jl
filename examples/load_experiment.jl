
# download the NMRData repository at https://github.com/AI4DBiological-Systems/NMRData

import NMRDataSetup

using FFTW
import PyPlot
import BSON

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])


### user inputs.

save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments"
project_name = "glucose-700"
experiment_full_path = "/home/roy/Documents/repo/NMRData/combination/BMRB-700/D-(+)-Glucose"
solvent_ppm_guess = 4.7
solvent_window_ppm = 0.1

### end inputs.

## load.
s_t, S, hz2ppmfunc, ppm2hzfunc, ν_0ppm, fs, SW, α_0ppm, β_0ppm, λ_0ppm, Ω_0ppm,
    results_0ppm, dic, α_solvent, β_solvent, λ_solvent, Ω_solvent,
    results_solvent = NMRDataSetup.loadspectrum(experiment_full_path;
    solvent_ppm = solvent_ppm_guess,
    solvent_window_ppm = solvent_window_ppm)

## store.
save_folder_path = joinpath(save_dir, project_name)
isdir(save_folder_path) || mkpath(save_folder_path); # make save folder if it doesn't exist.

save_path = joinpath(save_folder_path, "$(project_name).bson")
BSON.bson(save_path,
s_t = s_t,
fs = fs,
SW = SW,
ν_0ppm = ν_0ppm,
α_0ppm = α_0ppm,
β_0ppm = β_0ppm,
λ_0ppm = λ_0ppm)


## visualize.
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

offset_Hz = ν_0ppm - (ppm2hzfunc(0.3)-ppm2hzfunc(0.0))

N = length(s_t)
DFT_s = fft(s_t)
U_DFT, U_y, U_inds = NMRDataSetup.getwraparoundDFTfreqs(N, fs, offset_Hz)

# shift the DFT so that around 0 ppm (around `offset_Hz` Hz`) appears first in the frequency positions.
S_U = DFT_s[U_inds]
P_y = hz2ppmfunc.(U_y)



q = uu->NMRDataSetup.evalcomplexLorentzian(uu, α_0ppm, β_0ppm, λ_0ppm, 2*π*ν_0ppm)*fs
q_U = q.(U_y)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, real.(S_U), label = "data")
PyPlot.plot(P_y, real.(q_U), "--", label = "estimated 0ppm resonance component")
PyPlot.gca().invert_xaxis()

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("real part of data spectrum (shifted Discrete Fourier Transform)")


# It is difficult to fit resonance components perfectly in the frequency domain.
inds = findall(xx->(-0.1<xx<0.1), P_y)
P_y_near0 = P_y[inds]
S_U_near0 = S_U[inds]
q_U_near0 = q_U[inds]

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y_near0, real.(S_U_near0), label = "data: interpolated")
PyPlot.plot(P_y_near0, real.(q_U_near0), "--", label = "estimated 0ppm resonance component")
PyPlot.plot(P_y_near0, real.(S_U_near0), "x", label = "data: DFT points")

PyPlot.gca().invert_xaxis()

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("real part of data spectrum (shifted Discrete Fourier Transform)")


###normalized.

# # by model: normalize by the 0 ppm resonance component's estimate. sensitive to frequency mis-alignment.
# # Assumes the component has an α of 9, which corresponds to the 9 protons for DSS' 0ppm component.
# c = NMRDataSetup.evalcomplexLorentzian(ν_0ppm, 9.0, β_0ppm, λ_0ppm, 2*π*ν_0ppm)
# Z = q(ν_0ppm)/c

# by data.
# find closest U_y to ν_0ppm, then normalized according to S_U instead of q(blah freq)
val, ind = findmin( abs.(U_y .- ν_0ppm) )
Z = abs(S_U[ind])

y = S_U ./ Z

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, abs.(y))

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("abs")
PyPlot.title("magnitude of normalized data spectrum s.t. 0ppm peak has magnitude 1")
