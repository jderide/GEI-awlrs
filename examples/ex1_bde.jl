# Example 4.1: Brown, DeMarzo and Eaves economy.
# Two goods, three scenarios, agents A, B, B, two real contracts (the bond and a
# forward on good 2), no retention (A_s = 0), no transaction costs (D_0 = 0).
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using GEIawlrs, Random, LinearAlgebra

function build_bde(; eA = [2 0.5 1 1.5; 2 1 1 1.0], eB = [1 2.5 2 1.5; 1 2 2 2.0],
                   nB = 2, K = 5.7, strict = true)
    L, S = 2, 4
    pi = [1.0, 1/3, 1/3, 1/3]
    bA = repeat([0.25, 0.75], 1, S); bB = repeat([0.75, 0.25], 1, S)
    agents = vcat(Agent(bA, Float64.(eA), pi), [Agent(bB, Float64.(eB), pi) for _ in 1:nB])
    D = vcat([zeros(L, 2)], [Matrix(1.0I, L, 2) for _ in 2:S])
    A = [zeros(L, L) for _ in 1:S]
    Economy(agents, D, A; K = K, strict = strict)
end

function initial_prices_bde(rng)
    p = (0.6 .+ 0.4 .* rand(rng, 2, 4)) ./ 3
    p[1, 1] = 1.0; p[2, 1] = 0.6 + 0.2 * rand(rng)
    p
end

options_bde() = Options(gmax = 0.5, rho = 3.0, pmax0 = 1.5, pmax1 = 0.6, damping = 0.5)

if abspath(PROGRAM_FILE) == @__FILE__
    rec = solve_awlrs(build_bde(), initial_prices_bde(MersenneTwister(10)), options_bde())
    save_results(rec, joinpath(@__DIR__, "..", "results"), "ex1_bde")
end
