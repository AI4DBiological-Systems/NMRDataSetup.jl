function runNLopt!(  opt,
    p0::Vector{T},
    obj_func,
    grad_func,
    p_lbs,
    p_ubs;
    max_iters = 10000,
    xtol_rel = 1e-12,
    ftol_rel = 1e-12,
    maxtime = Inf) where T
#
    @assert length(p0) == length(p_lbs) == length(p_ubs)

    opt.maxeval = max_iters
    opt.lower_bounds = p_lbs
    opt.upper_bounds = p_ubs
    opt.xtol_rel = xtol_rel
    opt.ftol_rel = ftol_rel
    opt.maxtime = maxtime


    opt.min_objective = (xx, gg)->genericNLoptcostfunc!(xx, gg, obj_func, grad_func)

    # optimize.
    (minf, minx, ret) = NLopt.optimize(opt, p0)

    N_evals = opt.numevals

    return minf, minx, ret, N_evals
end

function genericNLoptcostfunc!(x::Vector{Float64}, df_x, f, df)::Float64

    #
    if length(df_x) > 0
    df_x[:] = df(x)
    end

    return f(x)
end