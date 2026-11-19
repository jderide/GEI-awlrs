# Example 4.2.1: Schmedders' variation of the endowments of agent A.
# The aggregate endowment of good 2 reaches 7 > K = 5.7, so the utility is not
# nondecreasing on the whole choice box (strict = false); the summary reports
# the largest utility index at the solution, which must stay below K.
include(joinpath(@__DIR__, "ex1_bde.jl"))

build_sch2() = build_bde(eA = [2 0.8 1.0 0.6; 2 2.4 3.0 1.8], strict = false)

if abspath(PROGRAM_FILE) == @__FILE__
    rec = solve_awlrs(build_sch2(), initial_prices_bde(MersenneTwister(10)), options_bde())
    save_results(rec, joinpath(@__DIR__, "..", "results"), "ex2_sch2")
end
