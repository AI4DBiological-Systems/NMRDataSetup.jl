
function evallinearitp(x0, y0, x1, y1, x::T)::T where T
    numerator = y0*(x1-x) + y1*(x-x0)
    return numerator/(x1-x0)
end


function scaletimeseries(s_data::Vector{Complex{T}}; new_max::T = convert(T, 10)) where T
    scale_factor = new_max/maximum( abs(s_data[n]) for n in eachindex(s_data) )
    s = collect( convert(Complex{T}, s_data[i]) * scale_factor for i in eachindex(s_data))

    return s, scale_factor
end
