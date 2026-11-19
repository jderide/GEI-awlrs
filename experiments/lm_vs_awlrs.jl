# Experiment: Levenberg-Marquardt started at the augmented Walrasian output (AW+LM)
# versus Levenberg-Marquardt started at the same initial prices (LM).
# Usage: julia --project=. experiments/lm_vs_awlrs.jl <tag> <first seed> <last seed>
# Writes experiments/results/<tag>_seeds<a>-<b>.csv
include(joinpath(@__DIR__, "..", "examples", "ex2_sch2.jl"))
include(joinpath(@__DIR__, "..", "examples", "ex2_sch3.jl"))
include(joinpath(@__DIR__, "..", "examples", "ex4.jl"))
using Printf
tag = ARGS[1]; s1 = parse(Int, ARGS[2]); s2 = parse(Int, ARGS[3])
# optional 4th argument: number of augmented Walrasian iterations (short warm start); then only that arm is run
KAW = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 0
quiet(o) = Options(; (f => getfield(o, f) for f in fieldnames(Options))..., verbose = false, polish = false)
# failed agent solves are retried once (Section 4 of the paper, Economy docstring)
withretry(ec) = (ec.retry = true; ec)
setup = Dict(
  "ex1_bde"  => (() -> withretry(build_bde()),  initial_prices_bde, quiet(options_bde())),
  "ex2_sch2" => (() -> withretry(build_sch2()), initial_prices_bde, quiet(options_bde())),
  "ex2_sch3" => (() -> withretry(build_sch3()), initial_prices_bde, quiet(options_bde())),
  "ex4"      => (() -> withretry(build_ex4(MersenneTwister(10))), initial_prices_ex4, quiet(options_ex4())))
build, start, opts = setup[tag]
if KAW > 0
  opts = Options(; (f => getfield(opts, f) for f in fieldnames(Options))..., maxit = KAW)
end
AWNAME = KAW > 0 ? "AW$(KAW)+LM" : "AW+LM"
const TOL = 1e-6
function report(io, tag, seed, m, k0, kaw, kf, naw, nlm, nfev, t, ec, pf, status)
  q = contract_prices(ec, pf); _, _, sol = excess(ec, pf); ec.nsolves -= length(sol)
  ba = any(any(s.zp .> ec.zbound - 1e-6) || any(s.zn .> ec.zbound - 1e-6) for s in sol)
  rk = rank(GEIawlrs.payoff_values(ec, pf); rtol = 1e-6)
  prel = join([@sprintf("%.5f", pf[2, n] / pf[1, n]) for n in 1:size(pf, 2)], ";")
  @printf(io, "%s,%d,%s,%.3e,%.3e,%.3e,%d,%d,%d,%d,%d,%.1f,%s,%d,%.5f,%.5f,%s,%s\n", tag, seed, m, k0,
          kaw, kf, kf <= TOL, naw, nlm, naw + nlm, nfev, t, ba, rk, q[1], q[2], prel, status)
end
mkpath(joinpath(@__DIR__, "results"))
file = joinpath(@__DIR__, "results", KAW > 0 ? @sprintf("%s_aw%d_seeds%d-%d.csv", tag, KAW, s1, s2) : @sprintf("%s_seeds%d-%d.csv", tag, s1, s2))
open(file, "w") do io
  println(io, "tag,seed,method,kappa0,kappa_aw,kappa_final,success,solves_aw,solves_lm,solves_total,lm_fev,time_s,bound_active,rank_M1,q0,q1,p_rel,status")
  for seed in s1:s2
    p0 = start(MersenneTwister(seed)); p0[1, 1] = 1.0
    # AW + LM
    ec = build()
    k0 = try
      ES, Z, _ = excess(ec, p0); certificate(p0, ES, Z)[1]
    catch
      NaN
    end
    if isnan(k0)   # an agent problem fails at the starting prices: start excluded for both methods
      for m in (AWNAME, "LM")
        KAW > 0 && m == "LM" && continue
        @printf(io, "%s,%d,%s,NaN,NaN,NaN,0,0,0,0,0,0.0,NaN,0,NaN,NaN,,start not evaluable\n", tag, seed, m)
      end
      flush(io); continue
    end
    ec.nsolves = 0
    t0 = time(); naw = 0
    try
      rec = solve_awlrs(ec, p0, opts); naw = ec.nsolves
      pol = refine(ec, rec[:awlrs_prices]; pmin = opts.pmin); nlm = ec.nsolves - naw
      pf = pol.cert < rec[:awlrs_certificate] ? pol.p : rec[:awlrs_prices]
      kf = min(pol.cert, rec[:awlrs_certificate])
      report(io, tag, seed, AWNAME, k0, rec[:awlrs_certificate], kf, naw, nlm, pol.nfev, time() - t0, ec, pf, "ok")
    catch err
      msg = replace(sprint(showerror, err), "," => ";", "\n" => " ")[1:min(end, 60)]
      n = ec.nsolves
      @printf(io, "%s,%d,%s,%.3e,NaN,NaN,0,%d,%d,%d,0,%.1f,NaN,0,NaN,NaN,,error: %s\n", tag, seed, AWNAME, k0,
              naw == 0 ? n : naw, naw == 0 ? 0 : n - naw, n, time() - t0, msg)
    end
    flush(io)
    KAW > 0 && continue
    # LM from the same start
    ec = build(); ec.nsolves = 0; t0 = time()
    try
      pol = refine(ec, p0; pmin = opts.pmin)
      report(io, tag, seed, "LM", k0, NaN, pol.cert, 0, ec.nsolves, pol.nfev, time() - t0, ec, pol.p, "ok")
    catch err
      msg = replace(sprint(showerror, err), "," => ";", "\n" => " ")[1:min(end, 60)]
      @printf(io, "%s,%d,LM,%.3e,NaN,NaN,0,0,%d,%d,0,%.1f,NaN,0,NaN,NaN,,error: %s\n", tag, seed, k0,
              ec.nsolves, ec.nsolves, time() - t0, msg)
    end
    flush(io)
  end
end
