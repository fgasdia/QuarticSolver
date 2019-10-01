"""
Finds roots of 4th order polynomials of the form
```math
ax⁴ + bx³ + cx² + dx + e = 0
```
"""
module QuarticSolver

"""
Return a cube roots of a complex number `z`.

See https://math.stackexchange.com/questions/394432/cube-roots-of-the-complex-numbers-1i
"""
function complexcbrts(z::Complex)
    r = abs(z)
    rroot = cbrt(r)
    θ = angle(z)
    dθ = θ/3

    r1 = rroot*cis(dθ)
    r2 = rroot*cis(2π/3 + dθ)
    r3 = rroot*cis(2π*2/3 + dθ)

    return r1, r2, r3
end

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

signedcbrt(x) = real(x) >= 0 ? cbrt(x) : -cbrt(-x)

"""
https://github.com/NKrvavica/fqs/blob/master/On_computing_roots.md
"""
function solvecubic(a0, b0, c0, d0)
    inva0 = 1/a0
    a, b, c = b0*inva0, c0*inva0, d0*inva0

    third = 1/3
    thirda = a*third
    thirda² = thirda^2
    sqrt3 = sqrt(3)

    f = third*b - thirda²
    g = thirda*(2*thirda² - b) + c
    h = 0.25*g^2 + f^3

    if f == g == h == 0
        r1 = -signedcbrt(c)
        return r1, r1, r1
    elseif h <= 0
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
        S = signedcbrt(-0.5*g + sqrth)
        U = signedcbrt(-0.5*g - sqrth)
        SplusU = S + U
        SminusU = S - U
        tmp = SminusU*sqrt3*0.5im
        r1 = SplusU - thirda
        r2 = -0.5*SplusU - thirda + tmp
        r3 = -0.5*SplusU - thirda - tmp
        return r1, r2, r3
    end
end

"""
https://github.com/NKrvavica/fqs
"""
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

    # s = sqrt(2*p + 2*real(z0) + 0im)
    s = sqrt(2*p + 2*z0)
    s == 0 ? t = z0^2 + r : t = -q/s

    r0, r1 = solvequadratic(1, s, z0 + t)
    r2, r3 = solvequadratic(1, -s, z0 - t)

    return r0-a₀, r1-a₀, r2-a₀, r3-a₀
end

function bookerquartic(B4, B3, B2, B1, B0)
    b3 = B3/(4*B4)
    b2 = B2/(6*B4)
    b1 = B1/(4*B4)
    b0 = B0/B4

    b3² = b3^2

    H = b2 - b3²
    Iv = b0 - 4*b3*b1 + 3*b2^2
    G = b1 - 3*b3*b2 + 2*b3²
    h = -Iv/12
    g = -G^2/4 - H*(H^2 + 3*h)
    tmpsqrt = sqrt(complex(g^2 + 4*h^3))
    p = (-g + tmpsqrt)/2
    p = ifelse(abs(p) > 1e-10, p, (-g - tmpsqrt)/2)

    s1, s2, s3 = complexcbrts(p)

    r1 = sqrt(s1 - h/s1 - H)
    r2 = sqrt(s2 - h/s2 - H)
    r3 = sqrt(s3 - h/s3 - H)
    r1 = ifelse(-2*r1*r2*r3/G ≈ 1, r1, -r1)

    q1 = r1 + r2 + r3 - b3
    q2 = r1 - r2 - r3 - b3
    q3 = -r1 + r2 - r3 - b3
    q4 = -r1 - r2 + r3 - b3

    return q1, q2, q3, q4
end

end # module
