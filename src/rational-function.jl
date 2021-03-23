
"""
    RationalFunction(p::AbstractPolynomial, q::AbstractPolynomial)
    p // q

Create a rational expression (`p/q`) from the two polynomials. 

There is no attempt to cancel common factors. The [`lowest_terms`](@ref) function attempts to do that.

For purposes of iteration, a rational function is treated like a tuple.

## Examples
```
julia> using Polynomials

julia> p,q = fromroots(Polynomial, [1,2,3]), fromroots(Polynomial, [2,3,4])
(Polynomial(-6 + 11*x - 6*x^2 + x^3), Polynomial(-24 + 26*x - 9*x^2 + x^3))

julia> pq = p // q
(-6 + 11*x - 6*x^2 + x^3) // (-24 + 26*x - 9*x^2 + x^3)

julia> lowest_terms(pq)
(-0.333333 + 0.333333*x) // (-1.33333 + 0.333333*x)

julia> pq(2.5)
-1.0

julia> pq(2) # uses first non-`0/0` ratio of `p⁽ᵏ⁾/q⁽ᵏ⁾`
-0.5

julia> pq^2
(36 - 132*x + 193*x^2 - 144*x^3 + 58*x^4 - 12*x^5 + x^6) // (576 - 1248*x + 1108*x^2 - 516*x^3 + 133*x^4 - 18*x^5 + x^6)

julia> derivative(pq)
(-108 + 180*x - 111*x^2 + 30*x^3 - 3*x^4) // (576 - 1248*x + 1108*x^2 - 516*x^3 + 133*x^4 - 18*x^5 + x^6)
```

!!! Note:
    The [RationalFunctions.jl](https://github.com/aytekinar/RationalFunctions.jl) was a helpful source of ideas.


"""
struct RationalFunction{T, X, P<:Polynomials.AbstractPolynomial{T,X}} <: AbstractRationalFunction{T,X,P}
    num::P
    den::P
    function RationalFunction(p::P, q::P) where {T,X, P<:Polynomials.AbstractPolynomial{T,X}}
        new{T,X,P}(p, q)
    end
    function RationalFunction(p::P, q::T) where {T,X, P<:Polynomials.AbstractPolynomial{T,X}}
        new{T,X,P}(p, q*one(P))
    end
    function RationalFunction(p::T, q::Q) where {T,X, Q<:Polynomials.AbstractPolynomial{T,X}}
        new{T,X,Q}(p*one(Q), q)
    end
end

RationalFunction(p,q)  = RationalFunction(promote(p,q)...)

# evaluation
(pq::RationalFunction)(x) = eval_rationalfunction(x, pq)

# Look like rational numbers
function Base.://(p::Polynomials.AbstractPolynomial,q::Polynomials.AbstractPolynomial)
    RationalFunction(p,q)
end


function Base.promote_rule(::Type{PQ}, ::Type{Q}) where {T,X,P,PQ<:RationalFunction{T,X,P},
                                                         S,Y,Q<:AbstractPolynomial{S,Y}}
    assert_same_variable(X,Y)
    R = promote_type(T,S)
    P′ = constructorof(P)
    RationalFunction{R,X,P′{R,X}}
end
function Base.promote_rule(::Type{PQ}, ::Type{PQ′}) where {T,X,P,PQ<:RationalFunction{T,X,P},
                                                           S,Y,Q,PQ′<:RationalFunction{S,Y,Q}}
    Polynomials.assert_same_variable(X,Y)
    R = promote_type(T,S)
    P′ = constructorof(promote_type(P,Q))
    RationalFunction{R,X,P′{R,X}}
end

function Base.promote_type(::Type{R}, ::Type{RR}) where {T,X,P,R<:RationalFunction{T,X,P},
                                                         S,  Q,RR<:RationalFunction{S,X,Q}}
    PQ = promote_type(P,Q)
    TS = promote_type(T,S)
    RationalFunction{TS, X, PQ}
end

function Base.convert(::Type{R}, pq::RR) where {T,X,P,R<:RationalFunction{T,X,P},
                                                S,  Q,RR<:RationalFunction{S,X,Q}}
    p,q = pq
    p′, q′ = convert(P, p), convert(P, q')
    RationalFunction(p′, q′)
end


function variable(::Type{PQ}) where {T,X,PQ <: RationalFunction{T,X}}
    x = variable(Polynomial)
    rational_function(PQ, x, one(x))
end
