
function extractsetting(
    ::Type{T},
    acq_lines::Vector{String},
    setting::String,
    ) where T

    query_string = "$(setting)="
    q = filter(xx->occursin(query_string, xx), acq_lines)

    if length(q) != 1
        println("Error: there isn't exactly one match for $setting in file.")
        return nothing
    end

    s = split(q[begin], "=")
    return parse(T, s[end])
end

# only for FID signals. Used for removing dead time at the start of the FID time series from Bruker spectrometers.
function estimatechangepoint(
    s::Vector{Complex{T}};
    threshold_factor::T = convert(T, 0.3), # trigger above 10 percent of max value.
    discard_factor::T = convert(T, 0.5), # discard 0.5*N extra samples
    verbose::Bool = true,
    ) where T <: AbstractFloat

    @assert zero(T) < discard_factor < one(T)
    @assert zero(T) < threshold_factor < one(T)

    max_s = maximum(abs(s[n]) for n in eachindex(s))
    threshold = threshold_factor*max_s

    N = 0

    loop_flag = true
    n = 0
    while n <= length(s) && loop_flag
        n += 1
        if abs(s[n]) > threshold
            loop_flag = false
            N = n
        end
    end
    if N == 0
        if verbose
            println("Warning: no dead time found. Returning original signal.")
        end
        return s, 0
    end

    k = round(Int, N + N*discard_factor)

    if k >= length(s)
        if verbose
            println("Warning: estimated truncation length is longer than the signal. Returning original signal.")
        end
        return s, 0
    end

    return s[k:end], k
end

abstract type NMRSettings end

struct Bruker1D1HSettings{T} <: NMRSettings
    TD::Int
    SW::T
    SFO1::T
    O1::T
    fs::T
end

# s is the offset/truncated and scaled version of s_data.
struct Data1D{T <: AbstractFloat, ST <: NMRSettings}
    s_data::Vector{Complex{T}}
    offset_ind::Int
    scale_factor::T
    s::Vector{Complex{T}}
    settings::ST
end

# front end.
function loadBruker(
    ::Type{T},
    load_path::String;
    FID_name::String = "fid",
    settings_file_name::String = "acqu",
    rescaled_max = convert(T, 10.0),
    ) where T <: AbstractFloat

    # acquisition parameters.
    acq_path = joinpath(load_path, settings_file_name)
    acq_lines = readlines(acq_path)

    TD = extractsetting(Int, acq_lines, "\$TD")
    SW = extractsetting(T, acq_lines, "\$SW")
    SFO1 = extractsetting(T, acq_lines, "\$SFO1")
    O1 = extractsetting(T, acq_lines, "\$O1")
    fs_acq = extractsetting(T, acq_lines, "\$SW_h")

    BYTORDA = extractsetting(T, acq_lines, "\$BYTORDA") # 1 for big-endian, 0 for small-endian.
    DTYPA = extractsetting(T, acq_lines, "\$DTYPA") # 0 for Int32, 2 for Float64.

    # get FID time series.
    binary_type = Int32
    if DTYPA == 2
        binary_type = Float64
    end

    FID_path = joinpath(load_path, FID_name)
    data_raw = reinterpret(Complex{binary_type}, read(FID_path))

    #endianness = :small
    data_binary = data_raw
    if BYTORDA == 1
        #endianness = :big
        data_binary .= ntoh.(data_raw)
    end

    s_data = convert(Vector{Complex{T}}, data_binary)

    #@show length(s_data)
    if length(s_data)*2 != TD

        N_acquired = round(Int, TD/2)
    
        test_segment = s_data[N_acquired+1:end]
        # sanity check.
        if !isempty(test_segment) && norm(test_segment) > eps(T)*10
            
            println("Error: length(s_data)*2 != TD, and the samples after TD/2 aren't zero.")
            println("The norm on those samples: ", norm(test_segment))
            @show length(s_data), length(s_data)*2, TD
            println("Truncate FID data at sample $N_acquired anyways.")
            println()
        end
        
        s_data = s_data[begin:N_acquired]
    end
    
    # z is the unscaled version of s.
    z, n0 = NMRDataSetup.estimatechangepoint(s_data; verbose = true)

    scale_factor =rescaled_max/maximum( abs(z[n]) for n in eachindex(z) )
    s = collect( convert(Complex{T}, z[i]) * scale_factor for i in eachindex(z))

    return Data1D(s_data, n0, scale_factor, s, Bruker1D1HSettings(TD, SW, SFO1, O1, fs_acq))
end