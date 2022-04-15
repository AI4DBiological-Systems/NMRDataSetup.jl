
# download the NMRData repository at https://github.com/AI4DBiological-Systems/NMRData

import NMRDataSetup

using FFTW
import PyPlot
import BSON

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

# creates the save path if it doesn't exist.
function loadbatchexperimentsinsubfolders(prefix_string::String, root_path::String, save_dir::String, solvent_ppm_guess::T, solvent_window_ppm::T) where T <: Real

    tmp = readdir(root_path, join = true)
    inds = findall(xx->isdir(xx), tmp)

    experiment_full_paths = tmp[inds]
    experiment_compound_names = readdir(root_path)[inds]
    project_names = collect( "$(prefix_string)-$(experiment_compound_names[i])" for i = 1:length(experiment_compound_names))

    isdir(save_dir) || mkdir(save_dir); # make save folder if it doesn't exist.

    for i = 1:length(experiment_compound_names)
        project_name = project_names[i]
        experiment_full_path = experiment_full_paths[i]

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

    end
end

### user inputs.
prefix_string = "BMRB-500-0.5mM" # the prefix string name that gets attached to each stored experiment.
root_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/similar_settings/BMRB-500-0.5mM/" # the folder that contain the experiment folders that you want to load.
save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments" # root path where the experiment BSON files will be stored.

solvent_ppm_guess = 4.7 # in units ppm.
solvent_window_ppm = 0.1 # in units ppm.
### end inputs.
loadbatchexperimentsinsubfolders(prefix_string, root_path, save_dir, solvent_ppm_guess, solvent_window_ppm)



### user inputs.
prefix_string = "BMRB-500-2mM" # the prefix string name that gets attached to each stored experiment.
root_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/similar_settings/BMRB-500-2mM/" # the folder that contain the experiment folders that you want to load.
save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments" # root path where the experiment BSON files will be stored.

solvent_ppm_guess = 4.7 # in units ppm.
solvent_window_ppm = 0.1 # in units ppm.
### end inputs.
loadbatchexperimentsinsubfolders(prefix_string, root_path, save_dir, solvent_ppm_guess, solvent_window_ppm)



### user inputs.
prefix_string = "BMRB-700-20mM" # the prefix string name that gets attached to each stored experiment.
root_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/BMRB/similar_settings/BMRB-700-20mM/" # the folder that contain the experiment folders that you want to load.
save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments" # root path where the experiment BSON files will be stored.

solvent_ppm_guess = 4.7 # in units ppm.
solvent_window_ppm = 0.1 # in units ppm.
### end inputs.
loadbatchexperimentsinsubfolders(prefix_string, root_path, save_dir, solvent_ppm_guess, solvent_window_ppm)



### user inputs.
prefix_string = "NRC-bioreactor-HEK293-2011" # the prefix string name that gets attached to each stored experiment.
root_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/NRC/bioreactor/HEK293/2011-06-06" # the folder that contain the experiment folders that you want to load.
save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments" # root path where the experiment BSON files will be stored.

solvent_ppm_guess = 4.7 # in units ppm.
solvent_window_ppm = 0.1 # in units ppm.
### end inputs.
loadbatchexperimentsinsubfolders(prefix_string, root_path, save_dir, solvent_ppm_guess, solvent_window_ppm)



### user inputs.
prefix_string = "NRC-bioreactor-2015" # the prefix string name that gets attached to each stored experiment.
root_path = "/home/roy/Documents/repo/NMRData/experiments_1D1H/NRC/bioreactor/Feng-2015-11-30" # the folder that contain the experiment folders that you want to load.
save_dir = "/home/roy/MEGAsync/outputs/NMR/experiments" # root path where the experiment BSON files will be stored.

solvent_ppm_guess = 4.7 # in units ppm.
solvent_window_ppm = 0.1 # in units ppm.
### end inputs.
loadbatchexperimentsinsubfolders(prefix_string, root_path, save_dir, solvent_ppm_guess, solvent_window_ppm)
