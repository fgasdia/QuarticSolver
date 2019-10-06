function sortroots(r1, r2, r3, r4)
    ar1, ar2, ar3, ar4 = abs(r1), abs(r2), abs(r3), abs(r4)
    ar1 > ar2 > ar3 > ar4 && return r1, r2, r3, r4
    ar1 > ar2 > ar4 > ar3 && return r1, r2, r4, r3
    ar1 > ar3 > ar2 > ar4 && return r1, r3, r2, r4
    ar1 > ar3 > ar4 > ar2 && return r1, r3, r4, r2
    ar1 > ar4 > ar2 > ar3 && return r1, r4, r2, r3
    ar1 > ar4 > ar3 > ar2 && return r1, r4, r3, r2
    ar2 > ar1 > ar3 > ar4 && return r2, r1, r3, r4
    ar2 > ar1 > ar4 > ar3 && return r2, r1, r4, r3
    ar2 > ar3 > ar1 > ar4 && return r2, r3, r1, r4
    ar2 > ar3 > ar4 > ar1 && return r2, r3, r4, r1
    ar2 > ar4 > ar1 > ar3 && return r2, r4, r1, r3
    ar2 > ar4 > ar3 > ar1 && return r2, r4, r3, r1
    ar3 > ar1 > ar2 > ar4 && return r3, r1, r2, r4
    ar3 > ar1 > ar4 > ar2 && return r3, r1, r4, r2
    ar3 > ar2 > ar1 > ar4 && return r3, r2, r1, r4
    ar3 > ar2 > ar4 > ar1 && return r3, r2, r4, r1
    ar3 > ar4 > ar1 > ar2 && return r3, r4, r1, r2
    ar3 > ar4 > ar2 > ar1 && return r3, r4, r2, r1
    ar4 > ar1 > ar2 > ar3 && return r4, r1, r2, r3
    ar4 > ar1 > ar3 > ar2 && return r4, r1, r3, r2
    ar4 > ar2 > ar1 > ar3 && return r4, r2, r1, r3
    ar4 > ar2 > ar3 > ar1 && return r4, r2, r3, r1
    ar4 > ar3 > ar1 > ar2 && return r4, r3, r1, r2
    ar4 > ar3 > ar2 > ar1 && return r4, r3, r2, r1
end

function refine(r1, r2, r3, r4)
    a, b, c, d = sortroots(r1, r2, r3, r4)

    # Initial estimates, 2 paths required
    α₀₁ = -real(r1 + r2)
    β₀₁ = real(r1*r2)

    α₀₂ = -real(r2 + r3)
    β₀₂ = real(r2*r3)

    γ₀₁, δ₀₁ = fastleastsquares(a, b, c, d, α₀₁, β₀₁)
    γ₀₂, δ₀₂ = fastleastsquares(a, b, c, d, α₀₂, β₀₂)

    α, β, γ, δ = optimize(a, b, c, d, α₀₁, β₀₁, γ₀₁, δ₀₁, α₀₂, β₀₂, γ₀₂, δ₀₂)
    r1, r2 = solvequadratic(1, α, β)
    r3, r4 = solvequadratic(1, γ, δ)

    return r1, r2, r3, r4
end

function optimize(a, b, c, d,
                  α₁, β₁, γ₁, δ₁,
                  α₂, β₂, γ₂, δ₂)
    numiter = 10
    ϵ₁ = totalerror(a, b, c, d)
    ϵ₂ = totalerror(a, b, c, d)
    ϵ₁vec = MVector{numiter,typeof(ϵ₁)}(undef)
    ϵ₂vec = MVector{numiter,typeof(ϵ₂)}(undef)

    e₁₁, e₁₂, e₁₃, e₁₄ = residual(a, b, c, d, α₁, β₁, γ₁, δ₁)
    e₂₁, e₂₂, e₂₃, e₂₄ = residual(a, b, c, d, α₂, β₂, γ₂, δ₂)

    iter = 1
    while iter <= numiter
        α₁, β₁, γ₁, δ₁ = backwardoptimize(a, b, c, d, α₁, β₁, γ₁, δ₁,
                                          e₁₁, e₁₂, e₁₃, e₁₄)
        α₂, β₂, γ₂, δ₂ = backwardoptimize(a, b, c, d, α₂, β₂, γ₂, δ₂,
                                          e₂₁, e₂₂, e₂₃, e₂₄)
        e₁₁, e₁₂, e₁₃, e₁₄ = residual(a, b, c, d, α₁, β₁, γ₁, δ₁)
        e₂₁, e₂₂, e₂₃, e₂₄ = residual(a, b, c, d, α₂, β₂, γ₂, δ₂)
        ϵ₁ = totalerror(e₁₁, e₁₂, e₁₃, e₁₄)
        ϵ₂ = totalerror(e₂₁, e₂₂, e₂₃, e₂₄)

        ϵ₁vec[iter] = ϵ₁
        ϵ₂vec[iter] = ϵ₂

        if iter > 4
            complete(ϵ₁, view(ϵ₁vec, iter-4:iter-1)) && return α₁, β₁, γ₁, δ₁
            complete(ϵ₂, view(ϵ₂vec, iter-4:iter-1)) && return α₂, β₂, γ₂, δ₂
        else
            ϵ₁ < 1e-15 && return α₁, β₁, γ₁, δ₁
            ϵ₂ < 1e-15 && return α₂, β₂, γ₂, δ₂
        end

        iter += 1
    end

    if ϵ₁ < ϵ₂
        return (α₁, β₁, γ₁, δ₁)
    else
        return (α₂, β₂, γ₂, δ₂)
    end
end

# complete(ϵ, ϵvec) = ϵ == 0 || any(ϵ == e for e in ϵvec)
complete(ϵ, ϵvec) = ϵ < 1e-15 || any(abs(ϵ - e) < 1e-15 for e in ϵvec)

function fastleastsquares(a, b, c, d, α₀, β₀)
    𝛷₁ = 1 + α₀^2 + β₀^2
    𝛷₂ = α₀*(1 + β₀)

    c₁ = a - α₀ + α₀*(b - β₀) + β₀*c
    c₂ = b - β₀ + α₀*c + β₀*d

    L₁ = sqrt(𝛷₁)
    L₃ = 𝛷₂/L₁
    L₂ = sqrt(𝛷₁ - 𝛷₂/𝛷₁*𝛷₂)

    y₁ = c₁/L₁
    y₂ = (c₂ - y₁*L₃)/L₂

    δ₀ = y₂/L₂
    γ₀ = (y₁ - δ₀*L₃)/L₁

    return γ₀, δ₀
end

totalerror(e₁, e₂, e₃, e₄) = abs(e₁) + abs(e₂) + abs(e₃) + abs(e₄)

function residual(a, b, c, d, α, β, γ, δ)
    e₁ = a - α - γ
    e₂ = b - β - α*γ - δ
    e₃ = c - β*γ - α*δ
    e₄ = d - β*δ

    return e₁, e₂, e₃, e₄
end

"""
Strobach 2010 The fast quartic solver
"""
function backwardoptimize(x₁, x₂, x₃, x₄, α, β, γ, δ,
                          e₁, e₂, e₃, e₄)
    U₂₃ = α - γ
    U₃₃ = β - δ - γ*U₂₃
    L₄₃ = -δ*U₂₃/U₃₃
    U₄₄ = β - δ - L₄₃*U₂₃

    x₁ = e₁
    x₂ = e₂ - γ*x₁
    x₃ = e₃ - δ*x₁ - γ*x₂
    x₄ = e₄ - δ*x₂ - L₄₃*x₃

    y₄ = x₄/U₄₄
    y₃ = (x₃ - U₂₃*y₄)/U₃₃
    y₂ = x₂ - U₂₃*y₃ - y₄
    y₁ = x₁ - y₃

    α += y₁
    β += y₂
    γ += y₃
    δ += y₄

    return α, β, γ, δ
end


computequartic(a, b, c, d, e, x) = a*x[1]^4 + b*x[2]^3 + c*x[3]^2 + d*x[4] + e
