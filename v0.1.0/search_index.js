var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "#MultivariateBases.AbstractPolynomialBasis",
    "page": "Introduction",
    "title": "MultivariateBases.AbstractPolynomialBasis",
    "category": "type",
    "text": "abstract type AbstractPolynomialBasis end\n\nPolynomial basis of a subspace of the polynomials [Section~3.1.5, BPT12].\n\n[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R. Semidefinite Optimization and Convex Algebraic Geometry. Society for Industrial and Applied Mathematics, 2012.\n\n\n\n\n\n"
},

{
    "location": "#MultivariateBases.FixedPolynomialBasis",
    "page": "Introduction",
    "title": "MultivariateBases.FixedPolynomialBasis",
    "category": "type",
    "text": "struct FixedPolynomialBasis{PT<:MP.AbstractPolynomialLike, PV<:AbstractVector{PT}} <: AbstractPolynomialBasis\n    polynomials::PV\nend\n\nPolynomial basis with the polynomials of the vector polynomials. For instance, FixedPolynomialBasis([1, x, 2x^2-1, 4x^3-3x]) is the Chebyshev polynomial basis for cubic polynomials in the variable x.\n\n\n\n\n\n"
},

{
    "location": "#MultivariateBases.MonomialBasis",
    "page": "Introduction",
    "title": "MultivariateBases.MonomialBasis",
    "category": "type",
    "text": "struct MonomialBasis{MT<:MP.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis\n    monomials::MV\nend\n\nMonomial basis with the monomials of the vector monomials. For instance, MonomialBasis([1, x, y, x^2, x*y, y^2]) is the monomial basis for the subspace of quadratic polynomials in the variables x, y.\n\n\n\n\n\n"
},

{
    "location": "#MultivariateBases.ScaledMonomialBasis",
    "page": "Introduction",
    "title": "MultivariateBases.ScaledMonomialBasis",
    "category": "type",
    "text": "struct ScaledMonomialBasis{MT<:MP.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis\n    monomials::MV\nend\n\nScaled monomial basis (see [Section 3.1.5, BPT12]) with the monomials of the vector monomials. Given a monomial x^alpha = x_1^alpha_1 cdots x_n^alpha_n of degree d = sum_i=1^n alpha_i, the corresponding polynomial of the basis is\n\nd choose alpha^frac12 x^alpha quad text where  quad\nd choose alpha = fracdalpha_1 alpha_2 cdots alpha_n\n\nFor instance, create a polynomial with the basis xy^2 xy creates the polynomial sqrt3 a xy^2 + sqrt2 b xy where a and b are new JuMP decision variables. Constraining the polynomial axy^2 + bxy to be zero with the scaled monomial basis constrains a/√3 and b/√2 to be zero.\n\n[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R. Semidefinite Optimization and Convex Algebraic Geometry. Society for Industrial and Applied Mathematics, 2012.\n\n\n\n\n\n"
},

{
    "location": "#MultivariateBases-1",
    "page": "Introduction",
    "title": "MultivariateBases",
    "category": "section",
    "text": "MultivariateBases.jl is a standardized API for multivariate polynomial bases based on the MultivariatePolynomials API.AbstractPolynomialBasis\nFixedPolynomialBasis\nMonomialBasis\nScaledMonomialBasis"
},

]}
