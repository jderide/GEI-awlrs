# Checks reported in Section 4, evaluated at the final prices of run_all.jl
# (last row of results/<tag>_history.csv), and the comparison of Table 4.
# Usage, after run_all.jl:  julia --project=. scripts/checks.jl
include(joinpath(@__DIR__, "..", "examples", "ex2_sch2.jl"))
include(joinpath(@__DIR__, "..", "examples", "ex2_sch3.jl"))
include(joinpath(@__DIR__, "..", "examples", "ex4.jl"))
using Printf

const RES = joinpath(@__DIR__, "..", "results")

"Final prices of a run: last row of the history file."
function final_prices(tag, L, S)
    lines = readlines(joinpath(RES, tag * "_history.csv"))
    head = split(lines[1], ","); last = parse.(Float64, split(lines[end], ","))
    p = zeros(L, S)
    for l in 1:L, s in 1:S
        p[l, s] = last[findfirst(==("p_$(l)_$(s)"), head)]
    end
    p
end

"kappa, rank of M1, largest c/m and w/m, largest position relative to M, largest |net position|."
function report(name, ec, p)
    ES, Z, sol = excess(ec, p)
    k, _ = certificate(p, ES, Z)
    rk = rank(GEIawlrs.payoff_values(ec, p); rtol = 1e-6)
    cm = maximum(maximum(s.c ./ ec.ub) for s in sol)
    wm = maximum(maximum(s.w ./ ec.ub) for s in sol)
    zM = maximum(max(maximum(s.zp), maximum(s.zn)) for s in sol) / ec.zbound
    zn = maximum(maximum(abs.(s.z)) for s in sol)
    @printf("%-6s kappa %.2e  rank M1 %d  max c/m %.3f  max w/m %.3f  max position/M %.3f  max |net position| %.3f\n",
            name, k, rk, cm, wm, zM, zn)
    return ES, Z, sol
end

println("Final prices of run_all.jl")
p = final_prices("ex1_bde", 2, 4); report("4.1", build_bde(), p)
p = final_prices("ex2_sch2", 2, 4); _, _, sol = report("4.2.1", build_sch2(), p)
ec = build_sch2()
idx = maximum(prod(s.c[:, n] .^ a.b[:, n]) for (s, a) in zip(sol, ec.agents) for n in 1:4)
@printf("       4.2.1: largest index prod_l c^alpha over agents and states %.2f (K = %.1f)\n", idx, ec.K)
p = final_prices("ex2_sch3", 2, 4); _, _, sol = report("4.2.2", build_sch3(), p)
@printf("       4.2.2: relative prices of good 2 in the second-stage states %s\n",
        join([@sprintf("%.3f", p[2, n] / p[1, n]) for n in 2:4], ", "))
@printf("       4.2.2: absolute contract positions between %.1f and %.1f\n",
        minimum(minimum(abs.(s.z)) for s in sol), maximum(maximum(abs.(s.z)) for s in sol))
rng = MersenneTwister(10); ec = build_ex4(rng)
p = final_prices("ex4", 6, 4); report("4.3", ec, p)
@printf("       4.3: a = (%s), b = (%s)\n",
        join([@sprintf("%.4f", a.e[1, 2] / a.e[1, 1]) for a in ec.agents], ", "),
        join([@sprintf("%.4f", a.e[1, 4] / a.e[1, 1]) for a in ec.agents], ", "))

# Table 4. Rows p_{s,1} (with p_{s,0} = 1), q by (7) and (8), interest rate, and
# the contract market imbalance max_j |Z_j| computed with this implementation.
# [14, Ex. 1] and [9, 5.1]: modified prices p~ built from the prices reported
# there, four digits. This paper: p and q from the final prices, and the
# imbalance at the final modified prices rounded to four digits.
println("\nTable 4 (Example 4.1)")
ec = build_bde()
pfin = final_prices("ex1_bde", 2, 4)
tab = [("[14, Ex. 1]", [1.0 0.2885 0.3057 0.3252; 0.7307 0.2224 0.2184 0.2138], nothing),
       ("[9, 5.1]",    [1.0 0.2890 0.3056 0.3256; 0.7375 0.2259 0.2193 0.2159], nothing),
       ("this paper",  pfin, round.(pfin, digits = 4))]
for (name, p, pz) in tab
    _, Z, _ = excess(ec, pz === nothing ? p : pz)
    q = contract_prices(ec, p)
    @printf("%-12s p1 = %s  q = %.4f, %.4f  rate %.2f%%  max|Z| %.2e\n", name,
            join([@sprintf("%.4f", p[2, n] / p[1, n]) for n in 1:4], " "), q[1], q[2],
            100 * (1 / q[1] - 1), maximum(abs.(Z)))
end
