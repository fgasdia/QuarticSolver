using Test
using Random
using Statistics
using BenchmarkTools
using PolynomialRoots

using QuarticSolver

#==
These tests simply compare PolynomialRoots using BigFloats to the QuarticSolver
using Float64 (and the Complex versions of both).
==#

# Generate random coefficients (PolynomialRoots requires coeffs in reverse order)
rng = MersenneTwister(235)

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

@testset "real, positive coefficient samples" begin
    pr_rresults = [sort(roots(collect(big.(coeffs))),
                        by=x->(round(real(x), digits=6), round(imag(x), digits=6))) for coeffs in eachcol(revrcoeffs)]
    qs_rresults = [sort([QuarticSolver.solvequartic(coeffs...)...],
                        by=x->(round(real(x), digits=6), round(imag(x), digits=6))) for coeffs in eachcol(rcoeffs)]

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
    maxabsdiffreal = maximum(maximum(abs.(x) for x in real.(rdiff)))
    maxabsdiffimag = maximum(maximum(abs.(x) for x in imag.(rdiff)))

    argmax(real.(rdiff))
    argmax(imag.(rdiff))
    argmin(real.(rdiff))
    argmin(imag.(rdiff))

    @test rfraction == 1
    @test all(avgdiffreal .< 1e-10)
    @test all(avgdiffimag .< 1e-10)
    @test maxabsdiffreal < 1e-8
    @test maxabsdiffimag < 1e-8
end

@testset "real, positive and negative coefficient samples" begin
    pr_rnresults = [sort(roots(collect(big.(coeffs))),
                         by=x->(round(real(x), digits=6), round(imag(x), digits=6))) for coeffs in eachcol(revrncoeffs)]
    qs_rnresults = [sort([QuarticSolver.solvequartic(coeffs...)...],
                         by=x->(round(real(x), digits=6), round(imag(x), digits=6))) for coeffs in eachcol(rncoeffs)]

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
    maxabsdiffreal = maximum(maximum(abs.(x) for x in real.(rndiff)))
    maxabsdiffimag = maximum(maximum(abs.(x) for x in imag.(rndiff)))

    argmax(imag.(rndiff))
    argmin(imag.(rndiff))

    @test rnfraction == 1
    @test all(avgdiffreal .< 1e-10)
    @test all(avgdiffimag .< 1e-10)
    @test maxabsdiffreal < 1e-8
    @test maxabsdiffimag < 1e-8
end

@testset "complex, positive coefficient samples" begin
    pr_cresults = [sort(roots(collect(big.(coeffs))),
                        by=x->(round(real(x), digits=6), round(imag(x), digits=6))) for coeffs in eachcol(revccoeffs)]
    qs_cresults = [sort([QuarticSolver.solvequartic(coeffs...)...],
                        by=x->(round(real(x), digits=6), round(imag(x), digits=6))) for coeffs in eachcol(ccoeffs)]

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
    maxabsdiffreal = maximum(maximum(abs.(x) for x in real.(cdiff)))
    maxabsdiffimag = maximum(maximum(abs.(x) for x in imag.(cdiff)))

    argmax(imag.(cdiff))
    argmin(imag.(cdiff))

    @test cfraction == 1
    @test all(avgdiffreal .< 1e-10)
    @test all(avgdiffimag .< 1e-10)
    @test maxabsdiffreal < 1e-8
    @test maxabsdiffimag < 1e-8
end

@testset "complex, positive and negative coefficient samples" begin
    pr_cnresults = [sort(roots(collect(big.(coeffs))),
                         by=x->(round(real(x), digits=6), round(imag(x), digits=6))) for coeffs in eachcol(revcncoeffs)]
    qs_cnresults = [sort([QuarticSolver.solvequartic(coeffs...)...],
                         by=x->(round(real(x), digits=6), round(imag(x), digits=6))) for coeffs in eachcol(cncoeffs)]

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
    maxabsdiffreal = maximum(maximum(abs.(x) for x in real.(cndiff)))
    maxabsdiffimag = maximum(maximum(abs.(x) for x in imag.(cndiff)))

    argmax(imag.(cndiff))
    argmin(imag.(cndiff))

    @test cnfraction == 1
    @test all(avgdiffreal .< 1e-10)
    @test all(avgdiffimag .< 1e-10)
    @test maxabsdiffreal < 1e-8
    @test maxabsdiffimag < 1e-8
end

@testset "complex, positive coefficients larger than 1" begin
    pr_lcresults = [sort(roots(collect(big.(coeffs))),
                         by=x->(round(real(x), digits=6), round(imag(x), digits=6))) for coeffs in eachcol(largerevccoeffs)]
    qs_lcresults = [sort([QuarticSolver.solvequartic(coeffs...)...],
                         by=x->(round(real(x), digits=6), round(imag(x), digits=6))) for coeffs in eachcol(largeccoeffs)]

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
    maxabsdiffreal = maximum(maximum(abs.(x) for x in real.(lcdiff)))
    maxabsdiffimag = maximum(maximum(abs.(x) for x in imag.(lcdiff)))

    argmax(imag.(lcdiff))
    argmin(imag.(lcdiff))

    @test lcfraction == 1
    @test all(avgdiffreal .< 1e-10)
    @test all(avgdiffimag .< 1e-10)
    @test maxabsdiffreal < 1e-8
    @test maxabsdiffimag < 1e-8
end

singleccoeffs = ccoeffs[:, 13]
singlerevccoeffs = reverse(singleccoeffs)
@benchmark roots($singlerevccoeffs)
@benchmark QuarticSolver.solvequartic($singleccoeffs[1],
                                      $singleccoeffs[2],
                                      $singleccoeffs[3],
                                      $singleccoeffs[4],
                                      $singleccoeffs[5])
