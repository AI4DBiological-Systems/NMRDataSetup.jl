



### load all experiments, then combine them using the first experiment's 0 ppm, sampling rate, and spectral width.
function loadexperiments(paths::Vector{String};
    weights = ones(Float64, length(paths)),
    solvent_ppm_guess = 4.7,
    solvent_window_ppm = 0.1)

    M = length(paths)
    s_t_set = Vector{Vector{Complex{Float64}}}(undef, M)
    ν_0ppm_set = Vector{Float64}(undef, M)
    fs_set = Vector{Float64}(undef, M)
    SW_set = Vector{Float64}(undef, M)

    for n = 1:M

        s_t, S, hz2ppmfunc, ppm2hzfunc, ν_0ppm_set[n], fs_set[n], SW_set[n], α_0ppm, β_0ppm, λ_0ppm, Ω_0ppm,
        results_0ppm, dic, α_solvent, β_solvent, λ_solvent, Ω_solvent,
        results_solvent = loadspectrum(paths[n];
        solvent_ppm = solvent_ppm_guess,
        solvent_window_ppm = solvent_window_ppm)
        
        s_t_set[n] = weights[n] .* s_t
    end

    return s_t_set, ν_0ppm_set, fs_set, SW_set
end

function assemblemixture(paths::Vector{String};
    weights = ones(Float64, length(paths)),
    solvent_ppm = 4.7,
    solvent_window_ppm = 0.1)

    s_t_set, ν_0ppm_set, fs_set, SW_set = loadexperiments(paths;
    weights = weights,
    solvent_ppm_guess = solvent_ppm,
    solvent_window_ppm = solvent_window_ppm)

    # combine via average.
    s_t = sum(s_t_set) ./ length(s_t_set)

    return s_t, ν_0ppm_set[1], fs_set[1], SW_set[1],
        s_t_set, ν_0ppm_set, fs_set, SW_set # just have it.
end



