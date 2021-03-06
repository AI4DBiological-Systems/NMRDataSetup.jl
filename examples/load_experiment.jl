
# download the NMRData repository at https://github.com/AI4DBiological-Systems/NMRData

import NMRDataSetup

using FFTW
import PyPlot
import BSON

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])


### user inputs.


solvent_ppm_guess = 4.7
solvent_window_ppm = 0.1

# save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments/NRC"
# project_name = "NRC-glucose-2018"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/NRC/misc/glucose/Sep-25-2018"

# save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments/NRC"
# project_name = "NRC-dmem-2012"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/NRC/misc/dmem_medium/Oct-22-2012"

# save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-600-100mM"
# project_name = "D-(+)-Glucose"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/glucose-600-bmse000855_1"

save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments/NRC"
project_name = "NRC-4_amino_acid-Jan2022-1"
experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/NRC/NRC_4_amino_acid_mixture_Jan_2022/1"

# save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments/misc"
# project_name = "bmse000297_ethanol"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/misc/bmse000297_ethanol"

# save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-500-100mM"
# project_name = "D-(+)-Glucose"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/similar_settings/BMRB-500-100mM/D-(+)-Glucose"
#
# save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments/misc"
# project_name = "bmse000038_glutamine"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/misc/bmse000038_glutamine"
#
# save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments/misc"
# project_name = "gissmo_900_leucine"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/misc/gissmo_900_leucine"
#
# save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments/misc"
# project_name = "gissmo_900_leucine_1"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/misc/gissmo_900_leucine/1"
#
# save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments/misc"
# project_name = "bmse000867_serine"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/misc/bmse000867_serine"
#
# save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments/BMRB-500-100mM"
# project_name = "L-Isoleucine"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/similar_settings/BMRB-500-100mM/L-Isoleucine"
#
# save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments/misc"
# project_name = "bmse000860_valine"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/misc/bmse000860_valine"
#
# save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments/misc"
# project_name = "bmse000915_methionine"
# experiment_full_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/misc/bmse000915_methionine"

# we're on tyrosine.

### end inputs.

## load.
isdir(save_dir) || mkdir(save_dir); # make save folder if it doesn't exist.

s_t, S, hz2ppmfunc, ppm2hzfunc, ??_0ppm, fs, SW, ??_0ppm, ??_0ppm, ??_0ppm, ??_0ppm,
    results_0ppm, dic, ??_solvent, ??_solvent, ??_solvent, ??_solvent,
    results_solvent = NMRDataSetup.loadspectrum(experiment_full_path;
    solvent_ppm = solvent_ppm_guess,
    solvent_window_ppm = solvent_window_ppm)

## store.
save_folder_path = joinpath(save_dir, project_name)
isdir(save_folder_path) || mkpath(save_folder_path); # make save folder if it doesn't exist.

save_path = joinpath(save_folder_path, "experiment.bson")
BSON.bson(save_path,
s_t = s_t,
fs = fs,
SW = SW,
??_0ppm = ??_0ppm,
??_0ppm = ??_0ppm,
??_0ppm = ??_0ppm,
??_0ppm = ??_0ppm,
dic = dic)


## visualize.
hz2ppmfunc = uu->(uu - ??_0ppm)*SW/fs
ppm2hzfunc = pp->(??_0ppm + pp*fs/SW)

offset_Hz = ??_0ppm - (ppm2hzfunc(0.3)-ppm2hzfunc(0.0))

N = length(s_t)
DFT_s = fft(s_t)
U_DFT, U_y, U_inds = NMRDataSetup.getwraparoundDFTfreqs(N, fs, offset_Hz)

# shift the DFT so that around 0 ppm (around `offset_Hz` Hz`) appears first in the frequency positions.
S_U = DFT_s[U_inds]
P_y = hz2ppmfunc.(U_y)



q = uu->NMRDataSetup.evalcomplexLorentzian(uu, ??_0ppm, ??_0ppm, ??_0ppm, 2*??*??_0ppm)*fs
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
# # Assumes the component has an ?? of 9, which corresponds to the 9 protons for DSS' 0ppm component.
# c = NMRDataSetup.evalcomplexLorentzian(??_0ppm, 9.0, ??_0ppm, ??_0ppm, 2*??*??_0ppm)
# Z = q(??_0ppm)/c

# by data.
# find closest U_y to ??_0ppm, then normalized according to S_U instead of q(blah freq)
val, ind = findmin( abs.(U_y .- ??_0ppm) )
Z = abs(S_U[ind])

y = S_U ./ Z

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_y, abs.(y))

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("abs")
PyPlot.title("magnitude of normalized data spectrum s.t. 0ppm peak has magnitude 1")
