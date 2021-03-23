using PolynomialRatios
using Polynomials
using Test

@testset "PolynomialRatios.jl" begin

    p,q = fromroots(Polynomial, [1,2,3]), fromroots(Polynomial, [2,3,4])
    r,s = ImmutablePolynomial([1,2,3], :x), ImmutablePolynomial([1,2,3], :y)
    t,u = Polynomial([1,2,pi]), ImmutablePolynomial([1.1, 2.2, 3.4], :y)
    
    # constructor
    @test p // q isa RationalFunction
    @test p // r isa RationalFunction
    @test_throws ArgumentError r // s
    
    pq = p // t # promotes to type of t
    @test pq isa RationalFunction{Float64, :x}
   
    # pieces
    pp,qq = p // q
    @test (pp,qq) == (p,q)
    
    # evaluation
    pq = p//q
    @test pq(2.5) ≈ p(2.5)/ q(2.5)
    @test pq(2) ≈ fromroots([1,3])(2) / fromroots([3,4])(2)
    
    # arithmetic
    rs = r // (r-1)
    x = 1.2
    @test (pq + rs)(x) ≈ (pq(x) + rs(x))
    @test (pq * rs)(x) ≈ (pq(x) * rs(x))
    @test (-pq)(x) ≈ -p(x)/q(x)
    
    # derivative
    pq = p // one(p)
    x = 2.3
    @test derivative(pq)(x) ≈ derivative(p)(x)
    pq = p // q
    dp,dq = derivative.((p,q))
    @test derivative(pq)(x) ≈ (dp(x)*q(x) - p(x)*dq(x))/q(x)^2
    
    # integral
    pq = p // (2*one(p))
    @test iszero(derivative(integrate(pq)) - pq)
    pq = p // q
    @test_throws ArgumentError integrate(pq)
    
    # lowest terms
    pq = p // q
    pp, qq = lowest_terms(pq)
    @test all(abs.(pp.(roots(qq))) .> 1/2)
end
