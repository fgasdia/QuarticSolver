using Test
using Random
using Statistics
using BenchmarkTools
using PolynomialRoots

using QuarticSolver
include("lwpc_quartic.jl")

rng = MersenneTwister(1234)

rcoeffs = rand(rng, Float64, 5, 1000)
rncoeffs = randn(rng, Float64, 5, 1000)
ccoeffs = rand(rng, ComplexF64, 5, 1000)
cncoeffs = randn(rng, ComplexF64, 5, 1000)
largeccoeffs = rand(rng, ComplexF64, 5, 1000)*100

revrcoeffs = reverse(rcoeffs, dims=1)
revrncoeffs = reverse(rncoeffs, dims=1)
revccoeffs = reverse(ccoeffs, dims=1)
revcncoeffs = reverse(cncoeffs, dims=1)
largerevccoeffs = reverse(largeccoeffs, dims=1)

# Real, uniform
pr_rresults = [sort(roots(collect(big.(coeffs))), by=x->(real(x), imag(x))) for coeffs in eachcol(revrcoeffs)]
qs_rresults = [sort([QuarticSolver.solvequartic(coeffs...)...], by=x->(real(x), imag(x))) for coeffs in eachcol(rcoeffs)]

rmatches = isapprox.(pr_rresults, qs_rresults, atol=1e-8)
rfraction = sum(rmatches)/length(rmatches)
rdiff = pr_rresults - qs_rresults
avgdiff = mean(rdiff)
avgdiffreal = mean(real.(rdiff))
avgdiffimag = mean(imag.(rdiff))
mindiffreal = minimum(real.(rdiff))
mindiffimag = minimum(imag.(rdiff))
maxdiffreal = maximum(real.(rdiff))
maxdiffimag = maximum(imag.(rdiff))

argmax(imag.(rdiff))
argmin(imag.(rdiff))

# Real, normal (includes negative)
pr_rnresults = [sort(roots(collect(big.(coeffs))), by=x->(real(x), imag(x))) for coeffs in eachcol(revrncoeffs)]
qs_rnresults = [sort([QuarticSolver.solvequartic(coeffs...)...], by=x->(real(x), imag(x))) for coeffs in eachcol(rncoeffs)]

rnmatches = isapprox.(pr_rnresults, qs_rnresults, atol=1e-8)
rnfraction = sum(rnmatches)/length(rnmatches)
rndiff = pr_rnresults - qs_rnresults
avgdiff = mean(rndiff)
avgdiffreal = mean(real.(rndiff))
avgdiffimag = mean(imag.(rndiff))
mindiffreal = minimum(real.(rndiff))
mindiffimag = minimum(imag.(rndiff))
maxdiffreal = maximum(real.(rndiff))
maxdiffimag = maximum(imag.(rndiff))

argmax(imag.(rndiff))
argmin(imag.(rndiff))

# Complex, uniform
pr_cresults = [sort(roots(collect(big.(coeffs))), by=x->(real(x), imag(x))) for coeffs in eachcol(revccoeffs)]
qs_cresults = [sort([QuarticSolver.solvequartic(coeffs...)...], by=x->(real(x), imag(x))) for coeffs in eachcol(ccoeffs)]

cmatches = isapprox.(pr_cresults, qs_cresults, atol=1e-8)
cfraction = sum(cmatches)/length(cmatches)
cdiff = pr_cresults - qs_cresults
avgdiff = mean(cdiff)
avgdiffreal = mean(real.(cdiff))
avgdiffimag = mean(imag.(cdiff))
mindiffreal = minimum(real.(cdiff))
mindiffimag = minimum(imag.(cdiff))
maxdiffreal = maximum(real.(cdiff))
maxdiffimag = maximum(imag.(cdiff))

argmax(imag.(cdiff))
argmin(imag.(cdiff))

# Complex, normal
pr_cnresults = [sort(roots(collect(big.(coeffs))), by=x->(real(x), imag(x))) for coeffs in eachcol(revcncoeffs)]
qs_cnresults = [sort([QuarticSolver.solvequartic(coeffs...)...], by=x->(real(x), imag(x))) for coeffs in eachcol(cncoeffs)]

cnmatches = isapprox.(pr_cnresults, qs_cnresults, atol=1e-8)
cnfraction = sum(cnmatches)/length(cnmatches)
cndiff = pr_cnresults - qs_cnresults
avgdiff = mean(cndiff)
avgdiffreal = mean(real.(cndiff))
avgdiffimag = mean(imag.(cndiff))
mindiffreal = minimum(real.(cndiff))
mindiffimag = minimum(imag.(cndiff))
maxdiffreal = maximum(real.(cndiff))
maxdiffimag = maximum(imag.(cndiff))

argmax(imag.(cndiff))
argmin(imag.(cndiff))

# Complex, larger than 1, uniform
pr_lcresults = [sort(roots(collect(big.(coeffs))), by=x->(real(x), imag(x))) for coeffs in eachcol(largerevccoeffs)]
qs_lcresults = [sort([QuarticSolver.solvequartic(coeffs...)...], by=x->(real(x), imag(x))) for coeffs in eachcol(largeccoeffs)]

lcmatches = isapprox.(pr_lcresults, qs_lcresults, atol=1e-8)
lcfraction = sum(lcmatches)/length(lcmatches)
lcdiff = pr_lcresults - qs_lcresults
avgdiff = mean(lcdiff)
avgdiffreal = mean(real.(lcdiff))
avgdiffimag = mean(imag.(lcdiff))
mindiffreal = minimum(real.(lcdiff))
mindiffimag = minimum(imag.(lcdiff))
maxdiffreal = maximum(real.(lcdiff))
maxdiffimag = maximum(imag.(lcdiff))

argmax(imag.(lcdiff))
argmin(imag.(lcdiff))

singleccoeffs = ccoeffs[:, 13]
singlerevccoeffs = reverse(singleccoeffs)
@benchmark roots($singlerevccoeffs)
@benchmark QuarticSolver.solvequartic($singleccoeffs[1],
                                      $singleccoeffs[2],
                                      $singleccoeffs[3],
                                      $singleccoeffs[4],
                                      $singleccoeffs[5])
