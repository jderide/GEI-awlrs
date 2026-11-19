# Runs the four examples of Section 4 and writes results/ (summary and history
# per example, and a table of all runs). Usage: julia --project=. run_all.jl
# Failed agent solves are not retried (Economy keyword retry = false, the default).
include(joinpath(@__DIR__, "examples", "ex2_sch2.jl"))
include(joinpath(@__DIR__, "examples", "ex2_sch3.jl"))
include(joinpath(@__DIR__, "examples", "ex4.jl"))
using Printf

const SEED = 10
out = joinpath(@__DIR__, "results")
rows = String[]
function run_example(tag, build, start, opts)
    rng = MersenneTwister(SEED)
    ec = build(rng)
    t = @elapsed rec = solve_awlrs(ec, start(rng), opts)
    save_results(rec, out, tag)
    push!(rows, @sprintf("%-9s %4d  %10.2e  %10.2e  %4d/%d  %-5s  %8d  %7.1f  %8.4f",
          tag, rec[:iterations], rec[:awlrs_certificate], rec[:certificate],
          rec[:rank_M1], length(rec[:q]), rec[:bound_active], rec[:nsolves], t,
          1 / rec[:q][1] - 1))
end
quiet(o) = Options(; (f => getfield(o, f) for f in fieldnames(Options))..., verbose = false)
run_example("ex1_bde", _ -> build_bde(), initial_prices_bde, quiet(options_bde()))
run_example("ex2_sch2", _ -> build_sch2(), initial_prices_bde, quiet(options_bde()))
run_example("ex2_sch3", _ -> build_sch3(), initial_prices_bde, quiet(options_bde()))
run_example("ex4", build_ex4, initial_prices_ex4, quiet(options_ex4()))
header = "example     it  cert(awlrs)  cert(final)  rank/J  bound  NLPsolves  time(s)  int.rate"
open(joinpath(out, "table.txt"), "w") do io
    println(io, header); foreach(r -> println(io, r), rows)
end
println(header); foreach(println, rows)
