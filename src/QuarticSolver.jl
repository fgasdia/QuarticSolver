#==
Finds roots of 4th order polynomials of the form
```math
ax⁴ + bx³ + cx² + dx + e = 0
```
based on the Fast Quartic and Cubic Solver
https://github.com/NKrvavica/fqs
by Nino Krvavica.
==#

__precompile__()

module QuarticSolver

# function solvequadratic(a0, b0, c0)
#     inva0 = 1/a0
#     a, b = b0*inva0, c0*inva0
#
#     a₀ = -a/2
#     Δ = a₀^2 - b
#     sqrtΔ = sqrt(complex(Δ))
#
#     r1 = a₀ - sqrtΔ
#     r2 = a₀ + sqrtΔ
#
#     return r1, r2
# end

"""
See Numerical Recipes 5.6, more numerically robust if A or C are small
"""
function solvequadratic2(a, b, c)
    Δ = sqrt(b^2 - 4*a*c)

    q = ifelse(real(conj(b)*Δ) >= 0, (b + Δ), (b - Δ))
    q /= -2

    if q == 0
        x₁ = x₂ = q
    else
        x₁ = q / a
        x₂ = c / q
    end

    return x₁, x₂
end


function solvecubic(a0, b0, c0, d0)
    inva0 = 1/a0
    a, b, c = b0*inva0, c0*inva0, d0*inva0

    third = 1//3
    thirda = a*third
    thirda² = thirda^2
    sqrt3 = sqrt(3)

    f = third*b - thirda²
    g = thirda*(2*thirda² - b) + c
    h = g^2/4 + f^3  # discriminant (`Δ` or `d` in some papers)

    if f == g == h == 0  # this check seems a little strong
        r1 = -cbrt(c)
        return r1, r1, r1
    elseif h isa Complex || h <= 0   # casus irreducibilis
        j = sqrt(-f)
        k = acos(-(g/2)/j^3)
        m = cos(third*k)
        n = sqrt3*sin(third*k)
        r1 = 2*j*m - thirda
        r2 = -j*(m + n) - thirda
        r3 = -j*(m - n) - thirda
        return r1, r2, r3
    else
        sqrth = sqrt(h)
        S = cbrt(-g/2 + sqrth)
        U = cbrt(-g/2 - sqrth)
        SplusU = S + U
        SminusU = S - U
        tmp = SminusU*sqrt3*im/2
        r1 = SplusU - thirda
        r2 = -SplusU/2 - thirda + tmp
        r3 = -SplusU/2 - thirda - tmp
        return r1, r2, r3
    end
end

"""
See numerical recipes 5.6.
Note we purposely exclude the direct solution for real a, b, c.
"""
function solvecubic2(a, b, c)
    a² = a^2
    thirda = a/3
    cbrti = im*cbrt(3)/2

    Q = (a² - 3b)/9
    R = (2*a²*a - 9*a*b + 27*c)/54

    Δ = sqrt(R^2 - Q^3)
    if real(conj(R)*Δ) >= 0
        A = -((R + Δ)^(1//3))
    else
        A = -((R - Δ)^(1//3))
    end

    B = ifelse(A == 0, zero(A), Q/A)

    x₁ = (A + B) - thirda
    x₂ = -(A + B)/2 - thirda + cbrti*(A - B)
    x₃ = -(A + B)/2 - thirda - cbrti*(A - B)

    return x₁, x₂, x₃
end

function solvequartic(a0, b0, c0, d0, e0)
    inva0 = 1/a0
    a, b, c, d = b0*inva0, c0*inva0, d0*inva0, e0*inva0

    a₀ = a/4
    a₀² = a₀^2

    # Subsidiary cubic equation
    p = 3*a₀² - b/2
    q = a*a₀² - b*a₀ + c/2
    r = 3*a₀²^2 - b*a₀² + c*a₀ - d

    # One root of cubic
    z0, _, _ = solvecubic(1, p, r, p*r-q^2/2)
    z0 = complex(z0)

    s = sqrt(2*p + 2*z0)
    s == 0 ? t = z0^2 + r : t = -q/s

    r0, r1 = solvequadratic(1, s, z0 + t)
    r2, r3 = solvequadratic(1, -s, z0 - t)

    return r0-a₀, r1-a₀, r2-a₀, r3-a₀
end

end # module
