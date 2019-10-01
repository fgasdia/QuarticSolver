"""
Finds roots of 4th order polynomials of the form
```math
ax⁴ + bx³ + cx² + dx + e = 0
```
"""
module QuarticSolver

using StaticArrays


"""
General closed-form solution of the quartic.

```math
ax⁴ + bx³ + cx² + dx + e = 0
```

Based on the ``general formula`` from https://en.wikipedia.org/wiki/Quartic_function
"""
function closed(a, b, c, d, e)

    a² = a^2
    b² = b^2

    p = (8*a*c - 3*b²)/(8*a²)
    q = (b^3 - 4*a*b*c + 8*a²*d)/(8*a^3)

    Δ₀ = c^2 - 3*b*d + 12*a*e
    Δ₁ = 2*c^3 - 9*b*c*d + 27*b²*e + 27*a*d^2 - 72*a*c*e

    Qroots = complexcbrts((Δ₁ + sqrt(Δ₁^2 - 4*Δ₀^3))/2)
    Qidx = 1  # try any of the roots

    Q = Qroots[Qidx]


    S = sqrt(-2/3*p + 1/(3*a)*(Q + Δ₀/Q))/2
    while isapprox(S, zero(S), atol=1e-9)
        # If S ≈ 0, choose a different cube root of Q
        Qidx += 1
        Q = Qroots[Qidx]
        S = sqrt(-2/3*p + 1/(3*a)*(Q + Δ₀/Q))/2
    end

    S² = S^2

    tmp = -b/(4*a)
    sqrt12 = sqrt(-4*S² - 2*p + q/S)/2
    sqrt34 = sqrt(-4*S² - 2*p - q/S)/2

    x1 = -tmp - S + sqrt12
    x2 = -tmp - S - sqrt12
    x3 = -tmp + S + sqrt34
    x3 = -tmp + S - sqrt34

    return x1, x2, x3, x4
end

function solvedepressedquartic(a, b, c, d)
    a3 = -b
    b3 = a*c - 4*d
    c3 = -a^2*d - c^2 + 4*b*d
end

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

"""
See https://www.mathemania.com/lesson/cardanos-formula-solving-cubic-equations/
"""
function solvecubic(a₃, a₂, a₁, a₀)
    # Depress
    a = a₂/a₃
    b = a₁/a₃
    c = a₀/a₃

    third = 1/3
    thirda = a*third
    sqrt3 = sqrt(3)

    p = b - a^2/3
    q = 2*a^3/27 - a*b/3 + c
    Δ = (q/2)^2 + (p/3)^3

    # Unit cube roots (in addition to `1`)
    ε₁ = -1/2 - im*sqrt3/2
    ε₂ = -1/2 + im*sqrt3/2

    if real(Δ) >= 0
        # 1 real, 2 complex roots
        v = cbrt(-q/2 + sqrt(Δ))
        w = -p/(3v)
        x1 = v - w
        x2 = v*ε₁ + w*ε₂
        x3 = v*ε₂ + w*ε₁
        return x1, x2, x3
    else  # Δ < 0
        # Complex roots
        r = abs(-q/2 + im*sqrt(-Δ))
        cbrtr = cbrt(r)
        ϕ = atan(-2*sqrt(-Δ)/q)
        ϕ < 0 && (ϕ = -ϕ)
        x1 = 2*cbrtr*cos(ϕ/3)
        x2 = 2*cbrtr*cos((ϕ+2π)/3)
        x3 = 2*cbrtr*cos((ϕ+4π)/3)
        return x1, x2, x3
    end
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

    cbrtr = cbrt(abs(p))
    thirdθ = angle(p)/3
    s1 = cbrtr*cis(thirdθ)
    s2 = cbrtr*cis(2π/3+thirdθ)
    s3 = cbrtr*cis(2π*2/3+thirdθ)

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
