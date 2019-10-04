"""
Finds roots of 4th order polynomials of the form
```math
ax⁴ + bx³ + cx² + dx + e = 0
```
based on the Fast Quartic and Cubic Solver
https://github.com/NKrvavica/fqs
by Nino Krvavica.
"""
module QuarticSolver

using StaticArrays

function solvequadratic(a0, b0, c0)
    inva0 = 1/a0
    a, b = b0*inva0, c0*inva0

    a₀ = -0.5*a
    Δ = a₀^2 - b
    sqrtΔ = sqrt(complex(Δ))

    r1 = a₀ - sqrtΔ
    r2 = a₀ + sqrtΔ

    return r1, r2
end

function solvecubic(a0, b0, c0, d0)
    inva0 = 1/a0
    a, b, c = b0*inva0, c0*inva0, d0*inva0

    third = 1/3
    thirda = a*third
    thirda² = thirda^2
    sqrt3 = sqrt(3)

    f = third*b - thirda²
    g = thirda*(2*thirda² - b) + c
    h = 0.25*g^2 + f^3  # discriminant (`Δ` or `d` in some papers)

    if f == g == h == 0  # this check seems a little strong
        r1 = -cbrt(c)
        return r1, r1, r1
    elseif h isa Complex || h <= 0   # casus irreducibilis
        j = sqrt(-f)
        k = acos(-0.5*g/j^3)
        m = cos(third*k)
        n = sqrt3*sin(third*k)
        r1 = 2*j*m - thirda
        r2 = -j*(m + n) - thirda
        r3 = -j*(m - n) - thirda
        return r1, r2, r3
    else
        sqrth = sqrt(h)
        S = cbrt(-0.5*g + sqrth)
        U = cbrt(-0.5*g - sqrth)
        SplusU = S + U
        SminusU = S - U
        tmp = SminusU*sqrt3*0.5im
        r1 = SplusU - thirda
        r2 = -0.5*SplusU - thirda + tmp
        r3 = -0.5*SplusU - thirda - tmp
        return r1, r2, r3
    end
end

function solvequartic(a0, b0, c0, d0, e0)
    inva0 = 1/a0
    a, b, c, d = b0*inva0, c0*inva0, d0*inva0, e0*inva0

    a₀ = 0.25*a
    a₀² = a₀^2

    # Subsidiary cubic equation
    p = 3*a₀² - 0.5*b
    q = a*a₀² - b*a₀ + 0.5*c
    r = 3*a₀²^2 - b*a₀² + c*a₀ - d

    # One root of cubic
    z0, _, _ = solvecubic(1, p, r, p*r-0.5*q^2)
    z0 = complex(z0)

    s = sqrt(2*p + 2*z0)
    s == 0 ? t = z0^2 + r : t = -q/s

    r0, r1 = solvequadratic(1, s, z0 + t)
    r2, r3 = solvequadratic(1, -s, z0 - t)

    return r0-a₀, r1-a₀, r2-a₀, r3-a₀
end

# function sortroots(r1, r2, r3, r4)
#     r2 > r1 && ((r1, r2) = (r2, r1))
#     r4 > r3 && ((r3, r4) = (r4, r3))
#     r3 > r2 && ((r2, r3) = (r3, r2))
#     return r1, r2, r3, r4
# end

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

    γ₀₁, δ₀₁ = fls(a, b, c, d, α₀₁, β₀₁)
    γ₀₂, δ₀₂ = fls(a, b, c, d, α₀₂, β₀₂)

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

    iter = 1
    while iter <= numiter
        α₁, β₁, γ₁, δ₁ = backwardoptimize(a, b, c, d, α₁, β₁, γ₁, δ₁)
        α₂, β₂, γ₂, δ₂ = backwardoptimize(a, b, c, d, α₂, β₂, γ₂, δ₂)

        e₁₁, e₁₂, e₁₃, e₁₄ = err(a, b, c, d, α₁, β₁, γ₁, δ₁)
        e₂₁, e₂₂, e₂₃, e₂₄ = err(a, b, c, d, α₂, β₂, γ₂, δ₂)
        ϵ₁ = totalerror(e₁₁, e₁₂, e₁₃, e₁₄)
        ϵ₂ = totalerror(e₂₁, e₂₂, e₂₃, e₂₄)

        ϵ₁vec[iter] = ϵ₁
        ϵ₂vec[iter] = ϵ₂

        if iter > 3
            complete(ϵ₁, view(ϵ₁vec, iter-3:iter)) && return α₁, β₁, γ₁, δ₁
            complete(ϵ₂, view(ϵ₂vec, iter-3:iter)) && return α₂, β₂, γ₂, δ₂
        else
            ϵ₁ == 0 && return α₁, β₁, γ₁, δ₁
            ϵ₂ == 0 && return α₂, β₂, γ₂, δ₂
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
complete(ϵ, ϵvec) = ϵ == 0 || ϵ in ϵvec

function fls(a, b, c, d, α₀, β₀)
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

function err(a, b, c, d, α, β, γ, δ)
    e₁ = a - α - γ
    e₂ = b - β - α*γ - δ
    e₃ = c - β*γ - α*δ
    e₄ = d - β*δ

    return e₁, e₂, e₃, e₄
end

"""
Strobach 2010 The fast quartic solver
"""
function backwardoptimize(x₁, x₂, x₃, x₄, α, β, γ, δ)
    e₁, e₂, e₃, e₄ = err(x₁, x₂, x₃, x₄, α, β, γ, δ)

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

end # module
