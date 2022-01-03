### load a few 700 MHz BMRB experiments, each with a known concentration.
# combine them into a single FID.


#include("NMRDataSetup.jl")
import NMRDataSetup

# using FFTW
# import PyPlot
# import BSON

# PyPlot.close("all")
# fig_num = 1

# PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

experiment_full_path = "/home/roy/Documents/data/NMR/NMRData/combination/BMRB-700/D-(+)-Glucose"
solvent_ppm_guess = 4.7
solvent_window_ppm = 0.1

s_t, S, hz2ppmfunc, ppm2hzfunc, ν_0ppm, fs, SW,
    α_0ppm, β_0ppm, λ_0ppm, Ω_0ppm,
       results_0ppm, dic, α_solvent, β_solvent, λ_solvent, Ω_solvent,
        results_solvent = NMRDataSetup.loadspectrum(experiment_full_path;
                                    solvent_ppm = solvent_ppm_guess,
                                    solvent_window_ppm = solvent_window_ppm)


@assert 1==4


# ### user inputs.
# BMRB_experiments_folder = "/home/roy/Documents/data/NMR/NMRData/combination"
# BMRB_setup_label = "BMRB-700-cal1"
# base_path_JLD = "/home/roy/Documents/data/NMR/NMRData/src/input/molecules"
# project_folder = "/home/roy/MEGAsync/outputs/NMR/combined/eight_700_cs_0.05"
# project_name = "eight_700"

# target_names = ["D-(+)-Glucose";
# "L-Histidine";
# "L-Isoleucine";
# "L-Phenylalanine";
# "L-Threonine";
# "L-Tryptophan";
# "L-Valine";
# "L-Serine"]
# ### end user inputs.

### exploration exp.
BMRB_experiments_folder = "/home/roy/Documents/data/NMR/NMRData/combination"
BMRB_setup_label = "BMRB-700"
base_path_JLD = "/home/roy/Documents/data/NMR/NMRData/src/input/molecules"
project_folder = "/home/roy/MEGAsync/outputs/NMR/combined/four_700_cs_0.05"
project_name = "four_700"

target_names = ["D-(+)-Glucose";
"L-Glutamine";
"L-Alanine";
"L-Serine"]
### exploration exp.

#target_names_path = "/home/roy/Documents/repo/NMRMetaboliteQuantifier/configs/DMEM_targeted_metabolites.txt"
#cs_config_path = "/home/roy/Documents/repo/NMRMetaboliteQuantifier/configs/cs_config.txt"
#base_path_GISSMO = "/home/roy/Documents/data/NMR/NMRData/src/input/GISSMO_data"
#project_name = "test_case"
#project_folder ="/home/roy/MEGAsync/outputs/NMR/dev"

### set up.
isdir(project_folder) || mkdir(project_folder); # make save folder if it doesn't exist.
BMRB_base_path = joinpath(BMRB_experiments_folder, BMRB_setup_label)
N_targets = length(target_names)
tmp = collect( NMRDataSetup.getGISSMOentry(target_names[i]) for i = 1:N_targets )
entries = collect( tmp[i].entry for i = 1:N_targets )
molecule_names = collect( tmp[i].molecule_name for i = 1:N_targets )

# get the list of experiment paths.
experiment_full_paths = collect( joinpath(BMRB_base_path, target_names[i]) for i = 1:length(target_names) )

@assert 5==4

### load all experiments, then combine them using the first experiment's 0 ppm, sampling rate, and spectral width.
function loadexperiments(experiment_full_paths, base_path_JLD)

    # load dummy compound, to speed things up.
    target_names = ["L-Serine";]
    N_targets = length(target_names)
    tmp = collect( NMRMetaboliteQuantifier.getGISSMOentry(target_names[i]) for i = 1:N_targets )
    entries = collect( tmp[i].entry for i = 1:N_targets )
    molecule_names = collect( tmp[i].molecule_name for i = 1:N_targets )

    M = length(experiment_full_paths)
    s_t_set = Vector{Vector{Complex{Float64}}}(undef, M)
    ν_0ppm_set = Vector{Float64}(undef, M)
    fs_set = Vector{Float64}(undef, M)
    SW_set = Vector{Float64}(undef, M)

    for n = 1:M
        s_t_set[n], _, _, _, ν_0ppm_set[n], fs_set[n], SW_set[n], _ = NMRMetaboliteQuantifier.runcalibration(experiment_full_paths[n],
            base_path_JLD,
            entries,
            molecule_names,
            NMRMetaboliteQuantifier.defaultconfig( true, false))
    end

    return s_t_set, ν_0ppm_set, fs_set, SW_set
end

s_t_set, ν_0ppm_set, fs_set, SW_set = loadexperiments(experiment_full_paths, base_path_JLD)

s_t = sum(s_t_set) ./ length(s_t_set)

####

_, _, hz2ppmfunc, ppm2hzfunc, ν_0ppm, fs, SW,
α_0ppm, λ_0ppm, Ω_0ppm, β_0ppm, data_dic,
as_set, Fs_set, updateαΩFSfunc_set,
p0_set, N_p_cs_sys_set, p0_group_set, as_set, Fs_set,
α_solvent, β_solvent, λ_solvent, Ω_solvent = NMRMetaboliteQuantifier.runcalibration(experiment_full_paths[1],
    base_path_JLD,
    entries,
    molecule_names,
    NMRMetaboliteQuantifier.defaultconfig( false, false))

s_t = s_t ./(maximum(abs.(s_t)))

offset_Hz = ν_0ppm - (ppm2hzfunc(0.3)-ppm2hzfunc(0.0))

N = length(s_t)
DFT_s = fft(s_t)
U_DFT, U, U_inds = NMRMetaboliteQuantifier.getwraparoundDFTfreqs(N, fs, offset_Hz)
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

save_path = joinpath(project_folder, "$(project_name).bson")
BSON.bson(save_path,
s_t = s_t,
p0_group_set = p0_group_set,
as_set = as_set,
Fs_set = Fs_set,
molecule_names = molecule_names,
fs = fs,
SW = SW,
ν_0ppm = ν_0ppm,
λ_0ppm = λ_0ppm,
β_0ppm = β_0ppm,
Ω_solvent = Ω_solvent,
λ_solvent = λ_solvent,
β_solvent = β_solvent)
