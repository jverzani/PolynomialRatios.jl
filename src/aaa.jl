## https://github.com/macd/BaryRational.jl/blob/main/src/aaa.jl
## https://arxiv.org/pdf/1612.00337.pdf for paper

import BaryRational: aaa, prz
export aaa, prz

"""
    BaryRationalFunction(xs, ys)
    fit(BaryRationalFunction, xs, ys)

Fit a rational function to the points `(xs,ys)` using the `aaa` method from [BaryRational](https://github.com/macd/BaryRational.jl/).

References: Nakatsukasa, Sete, and Trefethen; THE AAA ALGORITHM FOR RATIONAL APPROXIMATION [	arXiv:1612.00337](https://arxiv.org/pdf/1612.00337.pdf)

> Even on a disk or interval the algorithm may outperform existing methods, and on more complicated domains it is especially competitive. The core ideas are (1) representation of the rational approximant in barycentric form with interpolation at certain support points and (2) greedy selection of the support points to avoid exponential instabilities. The name AAA stands for "adaptive Antoulas--Anderson" in honor of the authors who introduced a scheme based on [Approximation of Large-Scale Dynamical Systems](https://doi.org/10.1137/1.9780898718713).

"""
struct BaryRationalFunction{T,X,P} <: AbstractRationalFunction{T,X,P}
    aaa::BaryRational.AAAapprox{Array{T,1}}
    function BaryRationalFunction(xs, ys; var=:x, kwargs...)
        o = BaryRational.aaa(collect(xs), collect(ys), kwargs...)
        T = eltype(o.x)
        X = var
        P = Polynomial{T, X}
        new{T,X,P}(o)
    end
end

## need a convert(BaryRationalFunction, p::RationalFunction) method...
BaryRationalFunction(p::AbstractPolynomial, q::AbstractPolynomial) =
    throw(ArgumentError("Method not defined for BaryRationalFunction"))

Polynomials.fit(::Type{<:BaryRationalFunction}, xs, ys; var=:x, kwargs...) =
    BaryRationalFunction(xs, ys; var=var, kwargs...)


Base.show(io::IO, F::BaryRationalFunction) = print(io, "BaryRational approximation")

Base.numerator(pq::BaryRationalFunction) = pqs(pq)[1]
Base.denominator(pq::BaryRationalFunction) = pqs(pq)[2]

# evaluation
(F::BaryRationalFunction)(x) = F.aaa(x)

function Base.convert(PQ::Type{RationalFunction}, F::BaryRationalFunction)
    P = eltype(PQ)
    z = variable(P) // one(P) 
    F.aaa(z)
end

function pqs(F::BaryRationalFunction{T,X,P}) where {T,X,P}
    pq = convert(RationalFunction, F)
    pqs(pq)
end








