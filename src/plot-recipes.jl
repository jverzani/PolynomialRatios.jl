## Plot recipe
## brute force this to avoid vertical asymptotes

@recipe function f(pq::AbstractRationalFunction{T}, a=nothing, b=nothing) where {T <: Real}

    p,q = lowest_terms(convert(RationalFunction,pq), method=:numerical)

    if a==nothing && b==nothing
        u=Polynomials.poly_interval(p)
        v=Polynomials.poly_interval(q)
        a,b = min(first(u), first(v)), max(last(u), last(v))
    end
    
    n = 251
    Δ = (b-a)/n

    # identify vertical asumptotes
    zs = filter(z -> a+Δ < z < b-Δ, sort(real.(filter(isreal, roots(q)))))
    xs = Any[]
    ys = Any[]
    if length(zs) > 0
        z₁ = popfirst!(zs)
        xs′ = range(a, stop=z₁ - Δ, length=50)
        ys′ = pq.(xs′)
        push!(xs, xs′)
        push!(ys, ys′)
        push!(xs, z₁)
        push!(ys, NaN)
        for z₂ ∈ zs
            xs′ = range(z₁ + Δ, stop=z₂ - Δ, length=50)
            ys′ = pq.(xs′)
            push!(xs, xs′)
            push!(ys, ys′)
            push!(xs, z₂)
            push!(ys, NaN)
            z₁ = z₂
        end
        xs′ = range(z₁+Δ, stop=b, length=50)
        ys′ = pq.(xs′)
        push!(xs, xs′)
        push!(ys, ys′)

        xs = vcat(xs...)
        ys = vcat(ys...)
    else
        xs = range(a, stop=b, length=251)
        ys = pq.(xs)
    end

    xs, ys
end

