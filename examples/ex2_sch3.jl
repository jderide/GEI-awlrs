# Example 4.2.2: Schmedders' economy with one agent of each type, so that the
# aggregate endowment is 3 for every good and state. At an equilibrium the
# relative prices of good 2 in the three scenarios cannot coincide (Section 4.2.2
# of the paper), but the computed prices are close to each other: M1(p) is close
# to losing rank and the contract positions are large. Neither the iterations nor
# the refinement reach an equilibrium; the run reports what they reach.
include(joinpath(@__DIR__, "ex1_bde.jl"))

build_sch3() = build_bde(nB = 1)

if abspath(PROGRAM_FILE) == @__FILE__
    rec = solve_awlrs(build_sch3(), initial_prices_bde(MersenneTwister(10)), options_bde())
    save_results(rec, joinpath(@__DIR__, "..", "results"), "ex2_sch3")
end
