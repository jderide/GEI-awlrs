# Example 4.3: heterogeneous agents, retention, and incomplete real contracts.
# Five agents, six goods, three scenarios (23 price unknowns).
# Contracts (J + 1 = 2 <= 3 scenarios, so M1(p) has full column rank generically):
#   contract 1: the bond, one unit of good 1 in every scenario;
#   contract 2: delivers good 2, in amounts (1.5, 0.5, 1.0) by scenario.
# No transaction costs (D_0 = 0).
# Retention: A_s = RET for every scenario; RET carries only the numeraire to the
# second stage. The utility has the retention term with weight rho_w and
# exponent gamma = 1/2, so it is concave and nondecreasing everywhere.
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using GEIawlrs, Random, LinearAlgebra

function retention_matrix()
    R = zeros(6, 6); R[1, 1] = 1.0
    R
end

function build_ex4(rng; RET = retention_matrix(), rho_w = 1e-2, gamma = 0.5)
    L, S = 6, 4
    pi = [1.0, 1/3, 1/3, 1/3]
    beta = [[3, 3, 2, 2, 1, 1], [3, 3, 2, 2, 1, 1], [2, 2, 3, 3, 1, 1],
            [2, 2, 3, 3, 1, 1], [1, 1, 1, 1, 6, 2]] ./ 12
    base = [[1, 4, 2, 2, 3, 1.5], [4, 1, 2, 2, 1, 1.5], [2, 1, 4, 3, 1, 1.5],
            [2, 4, 1, 3, 1, 1.5], [2, 2, 2, 2, 2, 2.0]]
    agents = Agent[]
    for i in 1:5
        e = repeat(base[i], 1, S)
        e[:, 2] .*= 1 + 0.2 * rand(rng)      # scenario 1: higher endowment
        e[:, 4] .*= 1 - 0.2 * rand(rng)      # scenario 3: lower endowment
        push!(agents, Agent(repeat(beta[i], 1, S), e, pi))
    end
    d2 = [1.5, 0.5, 1.0]                 # deliveries of contract 2 by scenario
    D = [zeros(L, 2) for _ in 1:S]
    for s in 2:S
        D[s][1, 1] = 1.0; D[s][2, 2] = d2[s - 1]
    end
    A = vcat([zeros(L, L)], [copy(RET) for _ in 2:S])
    Economy(agents, D, A; gamma = gamma, rho_w = rho_w)
end

# Start: all goods priced like the numeraire, state prices equal to the
# probabilities, with a 10% random perturbation.
function initial_prices_ex4(rng)
    p = (0.9 .+ 0.2 .* rand(rng, 6, 4)) ./ 3
    p[:, 1] .*= 3
    p[1, 1] = 1.0
    p
end

options_ex4() = Options(gmax = 1.0, pmax0 = 3.0, pmax1 = 1.5, r1 = 10.0, maxit = 30)

if abspath(PROGRAM_FILE) == @__FILE__
    rng = MersenneTwister(10)
    ec = build_ex4(rng)
    rec = solve_awlrs(ec, initial_prices_ex4(rng), options_ex4())
    save_results(rec, joinpath(@__DIR__, "..", "results"), "ex4")
end
