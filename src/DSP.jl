"""
    evalcomplexLorentzian(u, α, β::T, λ, Ω)::Complex{T} where T <: Real

Evaluates the complex Lorentzian function. Returns α*exp(im*β)/(λ + im*(2*π*u - Ω)).
"""
function evalcomplexLorentzian(u, α, β::T, λ, Ω)::Complex{T} where T <: Real
    return α*exp(im*β)/(λ + im*(2*π*u - Ω))
end

"""
    gettimerange(N::Int, fs::T) where T

Returns the time stamps for a sequence, starting at time 0. Returns zero(T):Ts:(N-1)*Ts, Ts = 1/fs.
"""
function gettimerange(N::Int, fs::T) where T
    Ts::T = 1/fs

    return zero(T):Ts:(N-1)*Ts
end

"""
    getwraparoundDFTfreqs(N::Int, fs::T, ν_begin::T) where T

Get the frequency indices for a sequence of length `N` that is sampled at `fs` Hz, with a starting frequency in `ν_begin` Hz.

...
# Inputs:


# Outputs, in the order returned.

- `U_DFT::LinRange{T,Int}`: frequency indices (in Hz) for `fft(x)`, where `length(x)` is `N` and values in `x` were sampled at `fs`.
- `U_y::Vector{T}`: frequency indices (in Hz) for the shited version of `fft(x)` that has the first frequency larger than `ν_begin` as the first entry. When fitting 1D 1H NMR data in the frequency domain, `U_y` is used instead of `U_DFT` with `ν_begin` being the 0 ppm frequency in Hz. One can then convert `U_y`, which is in Hz, to ppm, by using a conversion function appropriate for this experiment. See the `hz2ppmfunc` output from `loadspectrum()`.
- `U_inds::Vector{Int}`: `ff(x)[U_inds]` would give the correct spectrum values that correspond to the frequency indices in `U_y`.

Note that `U_y` doesn't contain the same numerical values as `U_DFT[U_inds]` because the latter is not neccessarily a monotonically increasing sequence, but `U_y` is. They both describe the same discrete-time frequencies though.
, , , , ,
    

# Examples
See load_experiment.jl in the /examples folder for an example.
"""
function getwraparoundDFTfreqs(N::Int, fs::T, ν_begin::T) where T

    U0 = getDFTfreqrange(N, fs)
    out, inds = wrapfreqrange(U0, ν_begin, fs)

    return U0, out, inds
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

function getDFTfreqrange(N::Int, fs::T)::LinRange{T} where T
    a = zero(T)
    b = fs-fs/N

    return LinRange(a, b, N)
end

  # case 1D.
function computeDTFTch3eq29(h, u::T, Λ)::Complex{T} where T <: Real

    # debug_array = zeros(Complex{T}, length(Λ))
    # debug_array2 = zeros(Complex{T}, length(Λ))
    # fs = 1/(Λ[2]-Λ[1])

    running_sum = zero(T)
    for i = 1:length(Λ)
        x = Λ[i]

        running_sum += h[i]*exp(-im*2*π*u*x)


        # # debug.
        # debug_array[i] = h[i]*exp(-im*2*π*u*x)
        #
        # #debug_array2[i] = h[i]*exp(-im*2*π*u*x +im*35*2*pi)
        # debug_array2[i] = h[i]*exp(-im*2*π*(u+3*fs)*x)
    end

    #println("debug: norm(debug_array -debug_array2 ) = ", norm(debug_array -debug_array2 ) )
    return running_sum
end
