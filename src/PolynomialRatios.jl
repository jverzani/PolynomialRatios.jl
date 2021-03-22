module PolynomialRatios

using LinearAlgebra, SparseArrays
using RecipesBase
import BaryRational

using Polynomials

export RationalFunction, BaryRationalFunction
export lowest_terms, poles, residues, partial_fraction

export pqs


include("common.jl")
include("rational-function.jl")

include("aaa.jl")
include("pade.jl")
include("plot-recipes.jl")

end
