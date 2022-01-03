
function loadBrukerfolderautophase(full_path::String)

    pdata_path = "$full_path/pdata/1"

    PyCall.py"""
    import numpy
    import nmrglue as ng

    s = $full_path
    procdata_dir = $pdata_path

    dic, data = ng.bruker.read(s)

    # ---time correction before FT
    precorr_time =  ng.bruker.remove_digital_filter(dic, data)
    precorr_frq0 = ng.proc_base.fft(precorr_time)
    precorr_frq = ng.proc_autophase.autops(precorr_frq0, 'acme')
    # precorr_frq = ng.proc_autophase.autops(precorr_frq0, 'peak_minima')

    #---time correction after FT
    postcorr_frq0 = ng.proc_base.fft(data)
    postcorr_frq1 = ng.bruker.remove_digital_filter(dic, postcorr_frq0, post_proc=True)
    postcorr_frq = ng.proc_autophase.autops(postcorr_frq1, 'acme')

    #---data processed using TopSpin (3.5pl7)
    bruker_frq = ng.bruker.read_pdata(procdata_dir)[1][::-1]

    guess_udic = ng.fileio.bruker.guess_udic(dic, data)
    """

    dic = PyCall.py"dic"
    data = PyCall.py"data"

    precorr_time = PyCall.py"precorr_time"
    precorr_frq0 = PyCall.py"precorr_frq0"
    precorr_frq = PyCall.py"precorr_frq"

    postcorr_frq0 = PyCall.py"postcorr_frq0"
    postcorr_frq1 = PyCall.py"postcorr_frq1"
    postcorr_frq = PyCall.py"postcorr_frq"

    bruker_frq = PyCall.py"bruker_frq"
    guess_udic = PyCall.py"guess_udic"

    return dic, data, precorr_time, precorr_frq0, precorr_frq,
            postcorr_frq0, postcorr_frq1, postcorr_frq, bruker_frq,
            guess_udic
end

function loadspectrum(full_path::String;
    solvent_ppm = 4.8,
    solvent_window_ppm = 0.3,
    N_Ω::Int = 20000,
    N_optim::Int = 1000,
    n_particles::Int = 3,
    max_iters::Int = 100,
    α_lower = 1e-6,
    α_upper_factor = 50.0,
    λ_lower = 1e-10,
    λ_upper = 20.0,
    λ_initial = 3.0,
    verbose_flag::Bool = false)

    dic, data, precorr_time, precorr_frq0, precorr_frq, postcorr_frq0, postcorr_frq1,
    postcorr_frq, Bruker_spectrums_nmrglue,
    guess_udic = loadBrukerfolderautophase(full_path)

    DE = dic["acqus"]["DE"] # in microseconds.

    BF1 = dic["acqus"]["BF1"]
    ## see Bruker TopSpin acquistion commands and parameters v 003.
    # TD - Time Domain; Number Of Raw Data Points. i.e., length(data)*2 - TD = 0.
    # total samples from both time-series.
    TD = dic["acqus"]["TD"]

    # SW - Spectral Width in ppm
    SW = dic["acqus"]["SW"]

    # SFO1 - SFO8 - Irradiation (carrier) Frequencies For Channels f1 to f8
    SFO1 = dic["acqus"]["SFO1"]

    O1 = dic["acqus"]["O1"]
    dic["acqus"]["NUC1"]

    fs_dic = dic["acqus"]["SW_h"]

    CAR = O1*1.0

    fs_dic = dic["acqus"]["SW_h"]

    N_diff = length(data)*2 - TD

    t_data = gettimerange(length(data), fs_dic)


    AQ = TD/(2*SW*SFO1)
    dead_time = t_data[end] - AQ

    st_ind0 = round(Int, dead_time/t_data[end] * length(data))

    st_ind_default = 100 # future: see if we can derive this from DE dead time later.
    st_ind = max(st_ind_default, st_ind0)

    s_t = data[st_ind:end]
    t = gettimerange(length(s_t), fs_dic)


    fs = (length(s_t)-1)/t[end] # this is really fs_dic.
    @assert abs(fs-fs_dic) < 1e-9


    α_0ppm, β_0ppm, λ_0ppm, Ω_0ppm,
    results_0ppm = estimatereferencecompoundfreq(CAR, fs, SW, s_t, t;
    N_Ω = N_Ω,
    N_optim = N_optim,
    n_particles = n_particles,
    max_iters = max_iters,
    α_lower = α_lower,
    α_upper_factor = α_upper_factor,
    λ_lower = λ_lower,
    λ_upper = λ_upper,
    λ_initial = λ_initial,
    verbose_flag = verbose_flag)

    ν_0ppm = Ω_0ppm ./ (2*π)

    hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
    ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

    # S is scaled DTFT of data.
    S = vv->(computeDTFTch3eq29(s_t, vv, t)/fs)

    # solvent.
    u_solvent_lower = ppm2hzfunc(solvent_ppm-solvent_window_ppm)
    u_solvent_upper = ppm2hzfunc(solvent_ppm+solvent_window_ppm)

    α_solvent, β_solvent, λ_solvent, Ω_solvent,
    results_solvent = estimatesinglet(fs, s_t, t,
                u_solvent_lower, u_solvent_upper;
                N_Ω = N_Ω,
                N_optim = N_optim,
                n_particles = n_particles,
                max_iters = max_iters,
                α_lower = α_lower,
                α_upper_factor = α_upper_factor,
                λ_lower = λ_lower,
                λ_upper = λ_upper,
                λ_initial = λ_initial,
                verbose_flag = verbose_flag)


    return s_t, S, hz2ppmfunc, ppm2hzfunc, ν_0ppm, fs, SW, α_0ppm, β_0ppm, λ_0ppm, Ω_0ppm,
    results_0ppm, dic, α_solvent, β_solvent, λ_solvent, Ω_solvent,
    results_solvent
end

function estimatereferencecompoundfreq( CAR::T,
    fs::T,
    SW::T,
    s_t::Vector{Complex{T}},
    t;
    N_Ω::Int = 20000,
    N_optim::Int = 1000,
    n_particles::Int = 3,
    max_iters::Int = 100,
    window_ppm::T = 0.3,
    α_lower::T = 1e-6,
    α_upper_factor::T = 50.0,
    λ_lower::T = 1e-10,
    λ_upper::T = 20.0,
    λ_initial::T = 3.0,
    verbose_flag::Bool = false) where T <: Real

    # prepare frequency range.
    ν0_initial, hz2ppmfunc0, ppm2hzfunc0 = getinitalguessreferencecompoundfreq(CAR, fs, SW)

    bp_a = ppm2hzfunc0(-window_ppm)
    bp_b = ppm2hzfunc0(window_ppm)

    return estimatesinglet(fs, s_t, t, bp_a, bp_b;
    N_Ω = N_Ω,
    N_optim = N_optim,
    n_particles = n_particles,
    max_iters = max_iters,
    α_lower = α_lower,
    α_upper_factor = α_upper_factor,
    λ_lower = λ_lower,
    λ_upper = λ_upper,
    λ_initial = λ_initial,
    verbose_flag = verbose_flag)

end

# λ_initial: m a value of 3 seems decent for the λ of DSS under 600 MHz.
# no error checking on whether λ_initial ∈ [λ_lower, λ_upper].
function estimatesinglet( fs::T,
    s_t::Vector{Complex{T}},
    t,
    u_lower,
    u_upper;
    N_Ω::Int = 20000,
    N_optim::Int = 1000,
    n_particles::Int = 3,
    max_iters::Int = 100,
    α_lower::T = 1e-6,
    α_upper_factor::T = 100.0,
    λ_lower::T = 1e-10,
    λ_upper::T = 20.0,
    λ_initial::T = 3.0,
    verbose_flag::Bool = false,
    xtol_rel::T = 1e-12,
    nlopt_seed::Int = 25) where T <: Real

    # prepare data.
    DTFT_s_scaled = vv->(computeDTFTch3eq29(s_t, vv, t)/fs)

    u_range_optim = LinRange(u_lower, u_upper, N_optim)
    S_U = DTFT_s_scaled.(u_range_optim)


    # set up optimization.
    max_val, max_ind = findmax(abs.(S_U))
    Ω = u_range_optim[max_ind] *2*π
    #println("Ω = ", Ω)
    costfunc = pp->evalcostreferencecompound(pp,
    S_U, u_range_optim, Ω)

    optim_lower_limit = [α_lower; -π; λ_lower]
    optim_upper_limit = [max_val*α_upper_factor; π; λ_upper]

    p0 = [max_val; 0.0; λ_initial]


    # optimize.
    op = Optim.Options( iterations = max_iters,
    store_trace = false,
    show_trace = verbose_flag)


    results = Optim.optimize( costfunc,
    p0,
    Optim.ParticleSwarm(; lower = optim_lower_limit,
        upper = optim_upper_limit,
        n_particles = n_particles),
    op)


    #α, β, λ = results.minimizer
    # println("α = ", α)
    # Q = [1221.2329692741823;
    # 1.5072497616902667;
    # 4.266537181116996]
    # println("cost Q = ", costfunc(Q))
    # println("cost p = ", costfunc(results.minimizer))

    # use NLopt.
    grad_func = xx->FiniteDiff.finite_difference_gradient(costfunc, xx)

    NLopt.srand(nlopt_seed)

    p_initial = results.minimizer
    #p_initial = p0

    opt = NLopt.Opt(:GN_ESCH, length(p_initial))
    #opt = NLopt.Opt(:LD_MMA, length(p_initial))

    min_cost_p, min_p, ret_p, numevals_p = runNLopt!(  opt,
        p_initial,
        costfunc,
        grad_func,
        optim_lower_limit,
        optim_upper_limit;
        max_iters = 100*max_iters,
        xtol_rel = xtol_rel)
    #
    opt = NLopt.Opt(:LD_MMA, length(min_p))
    min_cost_p, min_p, ret_p, numevals_p = runNLopt!(  opt,
        min_p,
        costfunc,
        grad_func,
        optim_lower_limit,
        optim_upper_limit;
        max_iters = 100*max_iters,
        xtol_rel = xtol_rel)

    α, β, λ = min_p

    # #
    # println("p_initial = ", p_initial)
    # println("α         = ", α)


    return α, β, λ, Ω, results
end

function evalcostreferencecompound( p::Vector{T},
    S_U::Vector{Complex{T}},
    u_range,
    Ω) where T <: Real

    # parse.
    α, β, λ = p

    cost = zero(T)
    for i = 1:length(u_range)
        q_u = evalcomplexLorentzian(u_range[i], α, β, λ, Ω)
        cost += abs2(q_u - S_U[i])
        #cost += (abs(q_u)-abs(S_U[i]))^2
    end

    return cost
end

function getinitalguessreferencecompoundfreq(CAR::T, fs::T, SW::T) where T <: Real

    ν0_initial = fs - CAR

    hz2ppmfunc0 = uu->(uu - ν0_initial)*SW/fs
    ppm2hzfunc0 = pp->(ν0_initial + pp*fs/SW)

    return ν0_initial, hz2ppmfunc0, ppm2hzfunc0
end

