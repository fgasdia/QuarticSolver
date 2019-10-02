# QuarticSolver
[![Build Status](https://travis-ci.com/EP-Guy/QuarticSolver.svg?branch=master)](https://travis-ci.com/EP-Guy/QuarticSolver) [![Build status](https://ci.appveyor.com/api/projects/status/0p7l4uanyyvnrxr8/branch/master?svg=true)](https://ci.appveyor.com/project/EP-Guy/quarticsolver/branch/master)

Find roots of 4th order polynomials (quartics) of the form:
```math
ax^4 + bx^3 + cx^2 + dx + e = 0
```

## Usage
```
root1, root2, root3 = solvequadratic(a, b, c, d, e)
```

Roots are unordered.

## References
This code is a Julia implementation of the Fast Quartic and Cubic Solver by Nino Krvavica (https://github.com/NKrvavica/fqs).
