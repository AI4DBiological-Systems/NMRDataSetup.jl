function evalcomplexLorentzian(u, α, β::T, λ, Ω)::Complex{T} where T <: Real
    return α*exp(im*β)/(λ + im*(2*π*u - Ω))
end

function gettimerange(N::Int, fs::T) where T
    Ts::T = 1/fs

    return zero(T):Ts:(N-1)*Ts
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