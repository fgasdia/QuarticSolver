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

function solvequadratic(a0, b0, c0)
    inva0 = 1/a0
    a, b = b0*inva0, c0*inva0

    a₀ = -0.5*a
    Δ = a₀^2 - b
    sqrtΔ = sqrt(Δ)

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

end # module
