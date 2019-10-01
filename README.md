# QuarticSolver
Find roots of 4th order polynomials (quartics) of the form:
```math
ax^4 + bx^3 + cx^2 + dx + e = 0
```

# Usage
```
root1, root2, root3 = solvequadratic(a, b, c, d, e)
```

Roots are unordered.

# References
This code is a Julia implementation of the Fast Quartic and Cubic Solver by Nino Krvavica (https://github.com/NKrvavica/fqs).
