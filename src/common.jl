import Polynomials: AbstractPolynomial, StandardBasisPolynomial
import Polynomials: isconstant, constantterm, assert_same_variable
import Polynomials: constructorof
import Polynomials: variable, indeterminate, derivative, integrate, ngcd

## ----
"""
    AbstractRationalFunction

Abstract type for holding ratios of polynomials.

Default methods for basic arithmetic operations are provided.

Numeric methods to cancel common factors, compute the poles, and return the residues are provided.
"""
abstract type AbstractRationalFunction{T,X,P} end


function Base.show(io::IO, pq::AbstractRationalFunction)
    p,q = pq
    print(io,"(")
    print(io,p)
    print(io, ") // (")
    print(io, q)
    print(io, ")")
end

## helper to make a rational function of given type
function rational_function(::Type{R}, p::P, q::Q) where {U,X, R<:AbstractRationalFunction{U,X},
                                                         T,   P<:AbstractPolynomial{T,X},
                                                         S,   Q<:AbstractPolynomial{S,X}}
    p′, q′ = promote(p,q)
    constructorof(R)(p,q)
end

Base.eltype(pq::Type{<:AbstractRationalFunction{T,X,P}}) where {T,X,P} = P
Base.eltype(pq::Type{<:AbstractRationalFunction{T,X}}) where {T,X} = Polynomial{T,X}
Base.eltype(pq::Type{<:AbstractRationalFunction{T}}) where {T} = Polynomial{T,:x}
Base.eltype(pq::Type{<:AbstractRationalFunction}) = Polynomial{Float64,:x}
Base.eltypeof(pq::AbstractRationalFunction{T,X,P}) where {T,X,P} = P

Base.size(F::AbstractRationalFunction) = ()
Base.isinf(F::AbstractRationalFunction) = false

## Look like rational numbers
# The p//q constructor is reserved for the `RationalFunction` type, but these should all be defined
# for other possible types.
function Base.://(p::PQ,q::PQ′) where {PQ <: AbstractRationalFunction, PQ′ <: AbstractRationalFunction}
    p0,p1 = p
    q0,q1 = p
    rational_function(promote_type(PQ, PQ′), p0*q1, p1*q0)
end

function Base.://(p::AbstractPolynomial,q::PQ) where {PQ <: AbstractRationalFunction}
    q0,q1 = q
    rational_function(PQ, p*q1, q0)
end
function Base.://(p::PQ, q::AbstractPolynomial) where {PQ <: AbstractRationalFunction}
    p0, p1 = p
    rational_function(PQ, p0, p1*q)
end
    
    
    
function Base.copy(pq::PQ) where {PQ <: AbstractRationalFunction}
    p,q = pq
    rational_function(PQ, p, q)
end


## ----

# Treat a RationalFunction as a tuple (num=p, den=q)
Base.length(pq::AbstractRationalFunction) = 2
function Base.iterate(pq, state=nothing)
    state == nothing && return (numerator(pq), 1)
    state == 1 && return (denominator(pq), 2)
    nothing
end
Base.collect(pq::AbstractRationalFunction{T,X,P}) where {T,X,P} = collect(P, pq)

# requires these field names
Base.numerator(pq::AbstractRationalFunction) = pq.num
Base.denominator(pq::AbstractRationalFunction) = pq.den

"""
    pqs(pq) 

Return `(p,q)`, where `pq=p/q`, as polynomials.
"""
pqs(pq::AbstractRationalFunction) = (numerator(pq), denominator(pq))

## ----

function Base.isnan(pq::AbstractRationalFunction) 
    p,q= pq
    iszero(p) && iszero(q)
end
function Base.iszero(pq::AbstractRationalFunction)
    p,q=pq
    iszero(p) && !iszero(q)
end
function Base.isone(pq::AbstractRationalFunction)
    p,q = pqs(pq)
    isconstant(p) && isconstant(q) && p == q
end

                                                             

## -----

indeterminate(pq::AbstractRationalFunction{T,X}) where {T,X} = X
isconstant(pq::AbstractRationalFunction; kwargs...) = all(Polynomials.isconstant.(lowest_terms(pq;kwargs...)))

function Base.zero(pq::R) where {R <: AbstractRationalFunction}
    p,q = pqs(pq)
    rational_function(R, zero(p), one(q))
end
function Base.one(pq::R) where {R <: AbstractRationalFunction}
    p,q = pqs(pq)
    rational_function(R, one(p), one(q))
end
function variable(pq::R) where {R <: AbstractRationalFunction}
    p,q = pqs(pq)
    rational_function(R, variable(p), one(q))
end

# use degree as largest degree of p,q after reduction
function Polynomials.degree(pq::AbstractRationalFunction)
    pq′ = lowest_terms(pq)
    maximum(degree.(pq′))
end

# Evaluation
function eval_rationalfunction(x, pq::AbstractRationalFunction)
    md = minimum(degree.(pq))
    num, den = pqs(pq)
    result = num(x)/den(x)
    while md >= 0
        !isnan(result) && return result
        num,den = derivative(num), derivative(den)
        result = num(x)/den(x)
        md -= 1
    end

    x*NaN
end
    
# Arithmetic
function Base.:-(pq::PQ) where {PQ <: AbstractRationalFunction}
    p, q = pqs(copy.(pq))
    rational_function(PQ, -p, q)
end

Base.:+(p::Number, q::AbstractRationalFunction) = q + p
Base.:+(p::AbstractRationalFunction,  q::Number) = p + q*one(p)
Base.:+(p::AbstractPolynomial, q::AbstractRationalFunction) = q + p
Base.:+(p::AbstractRationalFunction,  q::AbstractPolynomial) = p + (q//one(q))
Base.:+(p::AbstractRationalFunction, q::AbstractRationalFunction) = sum(promote(p,q))
function Base.:+(p::R, q::R) where {T,X,P,R <: AbstractRationalFunction{T,X,P}}
    p0,p1 = pqs(p)
    q0,q1 = pqs(q)
    rational_function(R, p0*q1 + p1*q0, p1*q1)
end

Base.:-(p::Number, q::AbstractRationalFunction) = -q +  p
Base.:-(p::AbstractRationalFunction,  q::Number) = p - q*one(p)
Base.:-(p::AbstractPolynomial, q::AbstractRationalFunction) = -q + p
Base.:-(p::PQ,  q::AbstractPolynomial) where {PQ <: AbstractRationalFunction} = p - rational_function(PQ,q, one(q))
function Base.:-(p::AbstractRationalFunction, q::AbstractRationalFunction)
    p′, q′ = promote(p,q)
    p′ - q′
end
function Base.:-(p::R, q::R) where {T,X,P,R <: AbstractRationalFunction{T,X,P}}
    p0,p1 = pqs(p)
    q0,q1 = pqs(q)
    rational_function(R, p0*q1 - p1*q0, p1*q1)
end

function Base.:*(p::Number, q::R) where {T, X, R <: AbstractRationalFunction{T,X}}
    q0,q1 = pqs(q)
    rational_function(R, (p*q0), q1)
end
function Base.:*(p::R,  q::Number) where {R <:AbstractRationalFunction}
    p0,p1 = pqs(p)
    rational_function(R, p0*q, p1)
end
function Base.:*(p::AbstractPolynomial, q::R) where {R <: AbstractRationalFunction}
    rational_function(R, p, one(p)) * q
end
Base.:*(p::R,  q::AbstractPolynomial) where {R <: AbstractRationalFunction} = p * rational_function(R,q, one(q))
Base.:*(p::AbstractRationalFunction, q::AbstractRationalFunction) = prod(promote(p,q))
function Base.:*(p::R, q::R) where {T,X,P,R <: AbstractRationalFunction{T,X,P}}
    p0,p1 = pqs(p)
    q0,q1 = pqs(q)
    rational_function(R, p0*q0, p1*q1)
end



function Base.:/(p::Number, q::R) where {R <: AbstractRationalFunction}
    q0,q1 = pqs(q)
    rational_function(R,p*q1, q0)
end
function Base.:/(p::R,  q::Number) where {R <: AbstractRationalFunction}
    p0,p1 = pqs(p)
    rational_function(R, p0, (p1*q))
end
Base.:/(p::AbstractPolynomial, q::PQ) where {PQ <: AbstractRationalFunction} = rational_function(PQ, p,one(p)) / q
function Base.:/(p::PQ,  q::AbstractPolynomial) where {PQ <: AbstractRationalFunction} 
    p0,p1 = pqs(p)
    rational_function(PQ,p0, p1*q)
end
function Base.:/(p::AbstractRationalFunction, q::AbstractRationalFunction)
    p′,q′ = promote(p,q)
    p/q
end
function Base.:/(p::PQ, q::PQ) where {T,X,P,PQ <: AbstractRationalFunction{T,X,P}}
    p0,p1 = pqs(p)
    q0,q1 = pqs(q)
    rational_function(PQ, p0*q1, p1*q0)
end

function Base.:^(pq::P, n::Int) where {P <: AbstractRationalFunction}
    p,q = pqs(pq)
    rational_function(P, p^n, q^n)
end

function Base.inv(pq::P) where {P <: AbstractRationalFunction}
    p,q = pqs(pq)
    rational_function(P, q, p)
end

# conj, transpose...

## derivative and integrals
function derivative(pq::P, n::Int=1) where {P <: AbstractRationalFunction}
    n <= 0 && return pq
    while n >= 1
        p,q = pqs(pq)
        pq = rational_function(P, (derivative(p)*q - p * derivative(q)), q^2)
        n -= 1
    end
    pq
end

function integrate(pq::P) where {P <: AbstractRationalFunction}
    p,q = pqs(pq)
    isconstant(q) && return rational_function(P, integrate(p), q)
    # XXX could work here, e.g.:
    # d,r = partial_fraction
    # ∫d + Σ∫r for each residue (issue with logs)
    throw(ArgumentError("Can only integrate rational functions with constant denominators"))
end

## ----
"""
    divrem(pq::AbstractRationalFunction; method=:numerical, kargs...)

Return `d,r` with `p/q = d + r/q` where `degree(numerator(r)) < degree(denominator(q))`, `d` a Polynomial, `r` a `AbstractRationalFunction`.

* `method`: passed to `gcd`
* `kwargs...`: passed to `gcd`
"""
function Base.divrem(pq::AbstractRationalFunction; method=:numerical, kwargs...)
    p,q = pqs(pq)
    degree(p) < degree(q) && return (zero(p), pq)

    d,r = divrem(p,q)
    (d, r//q)
    
end


# like Base.divgcd in rational.jl
# divide p,q by u
function _divgcd(v::Val{:euclidean}, pq; kwargs...)
    u = gcd(v, pqs(pq)...; kwargs...)
    p÷u, q÷u
end
function _divgcd(v::Val{:noda_sasaki}, pq; kwargs...)
    u = gcd(v, pqs(pq)...; kwargs...)
    p÷u, q÷u
end
function _divgcd(v::Val{:numerical}, pq; kwargs...)
    u,v,w,κ,θ = ngcd(pqs(pq)...; kwargs...) # u⋅v=p, u⋅w=q
    v, w
end


"""
    lowest_terms(pq::AbstractRationalFunctdion, method=:numerical)
    
Find GCD of `(p,q)`, `u`, and return `(p÷u)//(q÷u)`.
    
* `method`: passed to `gcd(p,q)`
* `kwargs`: passed to `gcd(p,q)`

"""
function lowest_terms(pq::PQ; method=:numerical, kwargs...) where {T,X,
                                                                   P<:StandardBasisPolynomial{T,X},
                                                                   PQ<:AbstractRationalFunction{T,X,P}}
    v,w = _divgcd(Val(method), pq; kwargs...)
    rational_function(PQ, v, w)
end

## ---- zeros, poles, ...
"""
    poles(pq::AbstractRationalFunction; method=:numerical, kwargs...)

For a rational function `p/q`, first reduces to lowest terms, then finds the roots and multiplicities of the resulting denominator.

"""
function poles(pq::AbstractRationalFunction; method=:numerical,  kwargs...)
    pq′ = lowest_terms(pq; method=method, kwargs...)
    den = denominator(pq′)
    mr = Polynomials.Multroot.multroot(den)
    (zs=mr.values, multiplicities = mr.multiplicities)
end

"""
    residues(pq::AbstractRationalFunction; method=:numerical,  kwargs...)

If `p/q =d + r/q`, returns `d` and the residues of a rational fraction `r/q`.

First expresses `p/q =d + r/q` with `r` of lower degree than `q` through `divrem`. 
Then finds the poles of `r/q`.
For a pole, `λj` of multiplicity `k` there are `k` residues, 
`rⱼ[k]/(z-λⱼ)^k`, `rⱼ[k-1]/(z-λⱼ)^(k-1)`, `rⱼ[k-2]/(z-λⱼ)^(k-2)`, …, `rⱼ[1]/(z-λⱼ)`.
The residues are found using this formula: 
`1/j! * dʲ/dsʲ (F(s)(s - λⱼ)^k` evaluated at `λⱼ` ([5-28](https://stanford.edu/~boyd/ee102/rational.pdf)).

There are several areas where numerical issues can arise. The `divrem`, the identification of multiple roots (`multroot`), the evaluation of the derivatives, ...

"""
function residues(pq::AbstractRationalFunction; method=:numerical,  kwargs...)

    
    d,r′ = divrem(pq)
    r = lowest_terms(r′; method=method, kwargs...)
    b,a = pqs(r)
    a′ = derivative(a)

    residues = Any[]
    mr = Polynomials.Multroot.multroot(a)

    for (λₖ, mₖ) ∈ zip(mr.values, mr.multiplicities)
        if mₖ == 1
            push!(residues, λₖ => [b(λₖ)/a′(λₖ)])
        else
            # returns rₖ,m, …, rₖ,1 where rₖ,i/(s-λₖ)ⁱ is part
            s = variable(a)
            F = lowest_terms(r*(s-λₖ)^mₖ)
            rs = [F(λₖ)]
            j! = 1
            for j ∈ 1:mₖ-1
                dF = lowest_terms(derivative(F))
                push!(rs, 1/j! * dF(λₖ))
                j! *= (j+1)
            end
            push!(residues, λₖ => rs)
        end
    end

    d, Dict(residues...)
end

"""
    partial_fraction(pq::AbstractRationalFunction; method=:numerical, kwargs...)

For a rational function `p/q = d + r/q`, with `degree(r) < degree(q)` returns `d` and
the terms that comprise `r/q`: For each pole with multiplicity returns
`rⱼ[k]/(z-λⱼ)^k`, `rⱼ[k-1]/(z-λⱼ)^(k-1)`, `rⱼ[k-2]/(z-λⱼ)^(k-2)`, …, `rⱼ[1]/(z-λⱼ)`.

Should be if `p/q` is in lowest terms and `d,r=partial_fraction(p//q)` that
`d + sum(r) - p//q ≈ 0`

"""
function partial_fraction(pq::AbstractRationalFunction; method=:numerical, kwargs...)
    d,r = residues(pq; method=method, kwargs...)
    s = variable(pq)
    d, partial_fraction(Val(:residue), r, s)
end

function partial_fraction(::Val{:residue}, r, s::PQ) where {PQ}
    terms = []
    for (λₖ,rsₖ) ∈ r
        for  (rⱼ,j) ∈ zip(rsₖ, length(rsₖ):-1:1)
            push!(terms, rⱼ/(s-λₖ)^j)
        end
    end
    terms

end
            
            
            
            
    

    





