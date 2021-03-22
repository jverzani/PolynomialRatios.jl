## ---- Fit
"""
    fit(::Type{RationalFunction}, r::Polynomial, m, n; var=:x)

Fit a Pade approximant ([`pade_fit`](@ref)) to `r`.

Examples:

```jldoctext
julia> using Polynomials, PolynomialRatios

julia> x = variable()
Polynomial(x)

julia> ex = 1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120 # Taylor polynomial for e^x
Polynomial(1.0 + 1.0*x + 0.5*x^2 + 0.16666666666666666*x^3 + 0.041666666666666664*x^4 + 0.008333333333333333*x^5)

julia> maximum(abs, exp(x) - fit(RationalFunction, ex, 1,1)(x) for x ∈ 0:.05:0.5)
0.017945395966538547

julia> maximum(abs, exp(x) - fit(RationalFunction, ex, 1,2)(x) for x ∈ 0:.05:0.5)
0.0016624471707165078

julia> maximum(abs, exp(x) - fit(RationalFunction, ex, 2,1)(x) for x ∈ 0:.05:0.5)
0.001278729299871717

julia> maximum(abs, exp(x) - fit(RationalFunction, ex, 2,2)(x) for x ∈ 0:.05:0.5)
7.262205147950951e-5
```
"""
function Polynomials.fit(::Type{RationalFunction},r::Polynomial, m::Integer, n::Integer;var=:x)
    p,q = pade_fit(r, m,n, var=var)
    p // q
end
    

## https://mathworld.wolfram.com/PadeApproximant.html
"""
    pade_fit(r::Polynomial, m,n)

For a polynomial `r` of degree `d ≥ m + n`, find a rational function `p/q` with
`degree(p) ≤ m`, `degree(q) ≤ n` and `q*r - p = x^{m+n+1}*s(x)` for some polynomial `s`.

This implementation sets up a system of equations to identify `p` and `q`.
"""
function pade_fit(p::Polynomial{T}, m::Integer, n::Integer; var=:x) where {T}
    d = degree(p)
    @assert (0 <= m) && (1 <= n) && (m + n <= d)

    # could be much more perfomant                
    c = convert(LaurentPolynomial, p) # for better indexing
    cs = [c[m+j-i] for j ∈ 1:n, i ∈ 0:n]
    
    qs′ = cs[:, 2:end] \ cs[:,1]
    qs = vcat(1, -qs′)

    cs = [c[0 + j - i] for j ∈ 0:m, i∈0:n]
    ps = cs * qs

    Polynomial(ps, var), Polynomial(qs,var)
end
export pade_fit

    
