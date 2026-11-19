"""
    GEIawlrs

Augmented Walrasian solver for two-stage exchange economies with incomplete
real financial markets and retention.

Companion code for "Solving equilibrium problems in economies with incomplete
financial markets and retention" (J. Deride, Set-Valued and Variational
Analysis).

Indexing: goods l = 1..L (good 1 is the numeraire), nodes s = 1..S (node 1 is
the first stage, nodes 2..S the second-stage scenarios), contracts j = 1..J
(contract 1 is the bond). Prices are the modified prices p~ of the paper, with
p[1,1] = 1 and contract prices q = sum_{s>1} D_s' p[:,s].
"""
module GEIawlrs

using JuMP, Ipopt, LinearAlgebra, Printf
import NLopt

export Agent, Economy, Options, excess, certificate, contract_prices,
       solve_awlrs, refine, save_results, bound_active

# ---------------------------------------------------------------- economy ---
"Agent data: Cobb-Douglas weights b[l,s], endowments e[l,s], node weights pi[s]."
struct Agent
    b::Matrix{Float64}
    e::Matrix{Float64}
    pi::Vector{Float64}
end

"""
    Economy(agents, D, A; K, rho_w, cmin, zbound)

`D[s]` is the L x J matrix of the contracts at node s: `D[1]` is the
transaction cost per unit sold (D_0 in the paper), `D[s]`, s > 1, the
deliveries. `A[s]`, s > 1, is the L x L retention matrix (`A[1]` is unused).
Utility of every agent, with x_s = prod_l c[l,s]^b[l,s] and
y_s = prod_l w[l,s]^b[l,s]:

    gamma = 0:      U(c, w) = - sum_s pi[s] (K - x_s - rho_w y_s)^2     (Examples 4.1, 4.2)
    gamma in (0,1): U(c, w) =   sum_s pi[s] (x_s^gamma + rho_w y_s^gamma)  (Example 4.3)

Choices are bounded by the box 0 <= c, w <= ub with
`ub[l,s] = 1.03 (aggregate endowment + retained aggregate endowment)`, so the
box is never active at a point where markets clear. For gamma = 0 the utility
is nondecreasing on the box only if `K >= (1 + rho_w) maximum(ub)`; the
constructor checks this (error, or a warning with `strict = false`; the
summary then reports the largest index x_s + rho_w y_s at the solution, to be
compared with K). For gamma in (0,1) the utility is concave and
nondecreasing everywhere.

`retry`: if true, an agent's problem on which Ipopt fails is solved once more
from the endowment and a zero portfolio (both solves are counted). If false
(default), the failure raises an error, which Phase II and the refinement turn
into a large penalty at that price. The runs of Table 3 (`run_all.jl`) use
`retry = false`, those of Section 4.4 (`experiments/`) `retry = true`.
"""
mutable struct Economy
    agents::Vector{Agent}
    D::Vector{Matrix{Float64}}
    A::Vector{Matrix{Float64}}
    K::Float64
    gamma::Float64
    rho_w::Float64
    cmin::Float64
    zbound::Float64
    ub::Matrix{Float64}
    retention::Bool
    models::Vector{Any}
    nsolves::Int
    retry::Bool
end

function Economy(agents, D, A; K = 5.7, gamma = 0.0, rho_w = 0.0, cmin = 1e-6,
                 zbound = 1e2, strict = true, retry = false)
    L, S = size(agents[1].e)
    etot = sum(a.e for a in agents)
    ub = copy(etot)
    for s in 2:S
        ub[:, s] .+= A[s] * etot[:, 1]
    end
    ub .*= 1.03
    retention = rho_w > 0 || any(!iszero, A[2:end])
    if gamma == 0 && K < (1 + rho_w) * maximum(ub)
        msg = "K = $K is below (1 + rho_w) max(ub) = $((1 + rho_w) * maximum(ub)): " *
              "the utility is not nondecreasing on the whole choice box"
        strict ? error(msg) : @warn(msg)
    end
    Economy(collect(agents), D, A, K, gamma, rho_w, cmin, zbound, ub, retention,
            Any[nothing for _ in agents], 0, retry)
end

ngoods(ec::Economy) = size(ec.ub, 1)
nnodes(ec::Economy) = size(ec.ub, 2)
ncontracts(ec::Economy) = size(ec.D[1], 2)

"No-arbitrage contract prices q = sum_{s>1} D_s' p[:,s] (equation (7) of the paper)."
contract_prices(ec::Economy, p) =
    sum(ec.D[s]' * p[:, s] for s in 2:nnodes(ec))

"Matrix M1(p) with rows p[:,s]' D_s, s > 1 (values of the contracts by scenario)."
payoff_values(ec::Economy, p) =
    reduce(vcat, [(p[:, s]' * ec.D[s]) for s in 2:nnodes(ec)])

# ----------------------------------------------------------- agent problem ---
function build_agent(ec::Economy, ag::Agent)
    L, S, J = ngoods(ec), nnodes(ec), ncontracts(ec)
    m = Model(Ipopt.Optimizer)
    set_silent(m)
    set_attribute(m, "tol", 1e-10)
    set_attribute(m, "acceptable_tol", 1e-8)
    set_attribute(m, "max_iter", 500)
    @variable(m, p[1:L, 1:S] in Parameter(1.0))
    @variable(m, q[1:J] in Parameter(0.0))
    @variable(m, ec.cmin <= c[l = 1:L, s = 1:S] <= ec.ub[l, s],
              start = clamp(ag.e[l, s], ec.cmin, ec.ub[l, s]))
    if ec.retention
        wlo = ec.rho_w > 0 ? ec.cmin : 0.0
        @variable(m, wlo <= w[l = 1:L, s = 1:S] <= ec.ub[l, s], start = max(wlo, 1e-2))
    else
        @variable(m, w[1:L, 1:S] == 0.0)
    end
    @variable(m, 0 <= zp[1:J] <= ec.zbound, start = 0.0)
    @variable(m, 0 <= zn[1:J] <= ec.zbound, start = 0.0)
    cd(s) = prod(c[l, s]^ag.b[l, s] for l in 1:L)
    wd(s) = prod(w[l, s]^ag.b[l, s] for l in 1:L)
    if ec.gamma > 0          # power utility, maximized
        if ec.rho_w > 0
            @objective(m, Max, sum(ag.pi[s] * (cd(s)^ec.gamma + ec.rho_w * wd(s)^ec.gamma) for s in 1:S))
        else
            @objective(m, Max, sum(ag.pi[s] * cd(s)^ec.gamma for s in 1:S))
        end
    else                     # quadratic utility with bliss level K
        if ec.rho_w > 0
            @objective(m, Min, sum(ag.pi[s] * (ec.K - cd(s) - ec.rho_w * wd(s))^2 for s in 1:S))
        else
            @objective(m, Min, sum(ag.pi[s] * (ec.K - cd(s))^2 for s in 1:S))
        end
    end
    # first-stage budget: transaction costs D_0 are charged on sales
    @constraint(m, sum(p[l, 1] * (c[l, 1] + w[l, 1] - ag.e[l, 1]
                    + sum(ec.D[1][l, j] * zn[j] for j in 1:J)) for l in 1:L)
                   + sum(q[j] * (zp[j] - zn[j]) for j in 1:J) <= 0)
    # second-stage budgets
    @constraint(m, [s = 2:S], sum(p[l, s] * (c[l, s] + w[l, s] - ag.e[l, s]
                    - sum(ec.D[s][l, j] * (zp[j] - zn[j]) for j in 1:J)
                    - sum(ec.A[s][l, k] * w[k, 1] for k in 1:L)) for l in 1:L) <= 0)
    return (model = m, p = p, q = q, c = c, w = w, zp = zp, zn = zn)
end

const _OK = (MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED, MOI.OPTIMAL)

"Solve the problem of agent i at prices p. Returns (c, w, zp, zn)."
function solve_agent!(ec::Economy, i::Int, p, q)
    if ec.models[i] === nothing
        ec.models[i] = build_agent(ec, ec.agents[i])
    end
    a = ec.models[i]
    set_parameter_value.(a.p, p)
    set_parameter_value.(a.q, q)
    optimize!(a.model)
    ec.nsolves += 1
    st = termination_status(a.model)
    if !(st in _OK) && ec.retry
        # retry once from a cold start: the warm start of the previous call can fail
        ag = ec.agents[i]
        set_start_value.(a.c, [clamp(ag.e[l, s], ec.cmin, ec.ub[l, s]) for l in axes(a.c, 1), s in axes(a.c, 2)])
        set_start_value.(a.zp, 0.0); set_start_value.(a.zn, 0.0)
        ec.retention && set_start_value.(a.w, ec.rho_w > 0 ? max(ec.cmin, 1e-2) : 1e-2)
        optimize!(a.model)
        ec.nsolves += 1
        st = termination_status(a.model)
    end
    st in _OK || error("agent $i: problem not solved ($st)")
    c, w, zp, zn = value.(a.c), value.(a.w), value.(a.zp), value.(a.zn)
    # warm start for the next call
    set_start_value.(a.c, c); set_start_value.(a.zp, zp); set_start_value.(a.zn, zn)
    ec.retention && set_start_value.(a.w, w)
    return c, w, zp, zn
end

# ------------------------------------------------------------ excess supply ---
"""
    excess(ec, p) -> (ES, Z, sol)

Excess supply `ES[l,s]` of the paper (individual excess supplies (2) summed
over agents, deliveries included) and aggregate net contract position
`Z = sum_i (z_i^+ - z_i^-)`.

Portfolio selection: each agent's net portfolio is replaced by the minimum-norm
portfolio with the same scenario values M1(p) z. It is the identity when M1(p)
has full column rank.
"""
function excess(ec::Economy, p)
    L, S, J = ngoods(ec), nnodes(ec), ncontracts(ec)
    q = contract_prices(ec, p)
    M1 = payoff_values(ec, p)
    P = pinv(M1; rtol = 1e-8) * M1
    ES = zeros(L, S); Z = zeros(J); sol = Any[]
    for i in eachindex(ec.agents)
        c, w, zp, zn = solve_agent!(ec, i, p, q)
        z = P * (zp - zn)
        ag = ec.agents[i]
        ES[:, 1] .+= ag.e[:, 1] .- c[:, 1] .- w[:, 1] .- ec.D[1] * zn
        for s in 2:S
            ES[:, s] .+= ag.e[:, s] .+ ec.D[s] * z .+ ec.A[s] * w[:, 1] .- c[:, s] .- w[:, s]
        end
        Z .+= z
        push!(sol, (c = c, w = w, z = z, zp = zp, zn = zn))
    end
    return ES, Z, sol
end

"""
    certificate(p, ES, Z) -> (value, (neg, comp, clr))

max( ||ES^-||_inf , max p .* |ES| , ||Z||_inf ): violation of ES >= 0, of
complementarity, and of contract clearing.
"""
function certificate(p, ES, Z)
    neg = maximum(max.(-ES, 0.0))
    comp = maximum(p .* abs.(ES))
    clr = isempty(Z) ? 0.0 : maximum(abs.(Z))
    return max(neg, comp, clr), (neg, comp, clr)
end

"True if some agent's portfolio is at the bound `zbound` at prices p."
function bound_active(ec::Economy, p; atol = 1e-6)
    _, _, sol = excess(ec, p)
    any(any(s.zp .> ec.zbound - atol) || any(s.zn .> ec.zbound - atol) for s in sol)
end

# ------------------------------------------------------- augmented Walrasian ---
"""
Algorithm options.

- `gmax`: bound of the dual domain G = [0, gmax]^(L S). All coordinates are
  included, the numeraire coordinate (1,1) too.
- `r0`, `r1`: first and last value of the proximal parameter r (geometric
  schedule). The augmented Walrasian is

      W_r(p, g) = min_{u in G} { <u, ES(p)> + |g - u|^2 / (2 r) } - rho |Z(p)|^2 ,

  with minimizer u = clip(g - r ES(p), 0, gmax). The parameter r of the paper
  (Section 4) is 1/r: r0 = 0.25 and r1 = 5 correspond to r^nu from 4 to 0.2.
  For every r > 0,
  inf_{g in G} W_r(p, g) = inf_{g in G} W(p, g), so W_r and W have the same
  maxinf-points.
- `rho`: weight of the contract-clearing penalty.
- `pmin`, `pmax0`, `pmax1`: price box (first-stage prices up to `pmax0`,
  second-stage modified prices up to `pmax1`).
- `tol`: stopping tolerance on the certificate; `maxit`: outer iterations.
- `ftol`, `xtol`, `maxeval`: Phase II (BOBYQA) relative tolerances and
  evaluation budget per price unknown.
- `damping`: step length of the price update.
- `polish`: refine the last iterate by a least-squares solve of (ES, Z) = 0.
"""
Base.@kwdef struct Options
    gmax::Float64 = 1.0
    r0::Float64 = 0.25
    r1::Float64 = 5.0
    rho::Float64 = 1.0
    pmin::Float64 = 1e-2
    pmax0::Float64 = 3.0
    pmax1::Float64 = 1.0
    tol::Float64 = 1e-2
    maxit::Int = 50
    ftol::Float64 = 1e-6
    xtol::Float64 = 1e-6
    maxeval::Int = 40
    damping::Float64 = 1.0
    polish::Bool = true
    verbose::Bool = true
end

pack(p) = vec(p)[2:end]
unpack(x, L, S) = reshape(vcat(1.0, x), L, S)

"Phase I: argmin over g in G of <g, ES>, i.e. g = gmax where ES < 0, else 0."
phase1(ES, gmax) = ifelse.(ES .< 0, gmax, 0.0)

"Value of the augmented Walrasian W_r(p, g) from ES(p), Z(p)."
function augmented_value(ES, Z, g, r, gmax, rho)
    u = clamp.(g .- r .* ES, 0.0, gmax)
    return sum(u .* ES .+ (g .- u) .^ 2 ./ (2r)) - rho * sum(abs2, Z)
end

"Phase II: maximize p -> W_r(p, g) over the price box by BOBYQA."
function phase2(ec::Economy, p, g, r, o::Options)
    L, S = ngoods(ec), nnodes(ec)
    n = L * S - 1
    lb = fill(o.pmin, n)
    ub = pack([s == 1 ? o.pmax0 : o.pmax1 for l in 1:L, s in 1:S])
    best = Ref(-Inf); bestx = Ref(clamp.(pack(p), lb, ub))
    function f(x, grad)
        v = try
            ES, Z, _ = excess(ec, unpack(x, L, S))
            augmented_value(ES, Z, g, r, o.gmax, o.rho)
        catch
            -1e6
        end
        if v > best[]
            best[] = v; bestx[] = copy(x)
        end
        return v
    end
    opt = NLopt.Opt(:LN_BOBYQA, n)
    NLopt.lower_bounds!(opt, lb); NLopt.upper_bounds!(opt, ub)
    NLopt.max_objective!(opt, f)
    NLopt.ftol_rel!(opt, o.ftol); NLopt.xtol_rel!(opt, o.xtol)
    NLopt.maxeval!(opt, o.maxeval * n)
    try
        NLopt.optimize(opt, clamp.(pack(p), lb, ub))
    catch
    end
    return unpack(bestx[], L, S)
end

"""
    solve_awlrs(ec, p0, options) -> record

Augmented Walrasian iterations from the prices `p0`, followed (if
`options.polish`) by the refinement step. The record holds the final prices,
contract prices, ES, Z, certificate, and the iteration history.
"""
function solve_awlrs(ec::Economy, p0, o::Options = Options())
    L, S = ngoods(ec), nnodes(ec)
    p = copy(p0); p[1, 1] = 1.0
    ES, Z, _ = excess(ec, p)
    cert, parts = certificate(p, ES, Z)
    hist = [(it = 0, cert = cert, neg = parts[1], comp = parts[2], clr = parts[3],
             r = o.r0, p = copy(p), ES = copy(ES), Z = copy(Z))]
    r = o.r0
    rg = (o.r1 / o.r0)^(1 / max(o.maxit - 10, 1))
    it = 0
    for k in 1:o.maxit
        cert <= o.tol && break
        g = phase1(ES, o.gmax)
        p1 = phase2(ec, p, g, r, o)
        p = o.damping .* p1 .+ (1 - o.damping) .* p
        ES, Z, _ = excess(ec, p)
        cert, parts = certificate(p, ES, Z)
        it = k
        push!(hist, (it = k, cert = cert, neg = parts[1], comp = parts[2],
                     clr = parts[3], r = r, p = copy(p), ES = copy(ES), Z = copy(Z)))
        o.verbose && @printf("it %3d  r=%.3g  cert=%.3e (ES-=%.2e comp=%.2e Z=%.2e)\n",
                             k, r, cert, parts...)
        r = min(r * rg, o.r1)
    end
    rec = Dict{Symbol,Any}(:iterations => it, :awlrs_certificate => cert,
                           :awlrs_prices => copy(p), :history => hist)
    if o.polish
        pol = refine(ec, p; pmin = o.pmin)
        o.verbose && @printf("refinement: certificate=%.2e  residual=%.2e  evaluations=%d\n",
                             pol.cert, pol.residual, pol.nfev)
        if pol.cert < cert
            p, ES, Z, cert = pol.p, pol.ES, pol.Z, pol.cert
            _, parts = certificate(p, ES, Z)
            push!(hist, (it = it + 1, cert = cert, neg = parts[1], comp = parts[2],
                         clr = parts[3], r = r, p = copy(p), ES = copy(ES), Z = copy(Z)))
        end
        rec[:refinement] = pol
    end
    _, _, sol = excess(ec, p)
    merge!(rec, Dict(:prices => p, :q => contract_prices(ec, p), :ES => ES, :Z => Z,
                     :certificate => cert, :allocations => sol,
                     :bound_active => any(any(s.zp .> ec.zbound - 1e-6) ||
                                          any(s.zn .> ec.zbound - 1e-6) for s in sol),
                     :rank_M1 => rank(payoff_values(ec, p); rtol = 1e-6),
                     :nsolves => ec.nsolves))
    if ec.gamma == 0
        idx = maximum(prod(s.c[:, n] .^ ec.agents[i].b[:, n]) +
                      ec.rho_w * prod(s.w[:, n] .^ ec.agents[i].b[:, n])
                      for (i, s) in enumerate(sol), n in 1:S)
        rec[:bliss] = (idx, ec.K)
    end
    return rec
end

# ---------------------------------------------------------------- refinement ---
"""
    refine(ec, p0; tol, maxit, pmin)

Levenberg-Marquardt iterations on F(p) = (ES(p), Z(p)) = 0 with a
forward-difference Jacobian, started at `p0`. It is a local method for interior
equilibria (positive prices, ES = 0) where ES and Z are single-valued and
smooth; it is not covered by the convergence theory of the paper.
"""
function refine(ec::Economy, p0; tol = 1e-9, maxit = 300, pmin = 1e-4, h = 1e-5)
    L, S = ngoods(ec), nnodes(ec)
    nfev = Ref(0)
    function F(x)
        nfev[] += 1
        try
            ES, Z, _ = excess(ec, unpack(x, L, S))
            return vcat(vec(ES), Z)
        catch
            return fill(1e3, L * S + ncontracts(ec))
        end
    end
    x = max.(pack(p0), pmin)
    f = F(x); lam = 1e-3
    Jm = zeros(length(f), length(x)); fresh = false
    for _ in 1:maxit
        maximum(abs, f) <= tol && break
        if !fresh
            for k in eachindex(x)
                xk = copy(x); dk = h * max(abs(x[k]), 1e-2); xk[k] += dk
                Jm[:, k] = (F(xk) .- f) ./ dk
            end
            fresh = true
        end
        H = Jm' * Jm
        d = -(H + lam * Diagonal(diag(H) .+ 1e-12)) \ (Jm' * f)
        xn = max.(x .+ d, pmin)
        fn = F(xn)
        if norm(fn) < norm(f)
            small = norm(xn - x) <= 1e-13 * (1 + norm(x))
            x, f, lam, fresh = xn, fn, max(lam / 10, 1e-12), false
            small && break
        else
            lam *= 10
            lam > 1e10 && break
        end
    end
    p = unpack(x, L, S)
    ES, Z, _ = excess(ec, p)
    cert, _ = certificate(p, ES, Z)
    return (p = p, ES = ES, Z = Z, cert = cert, residual = norm(f), nfev = nfev[])
end

# ------------------------------------------------------------------- output ---
fmt(M::AbstractMatrix) = join([join([@sprintf("%10.5f", v) for v in M[i, :]], " ")
                               for i in axes(M, 1)], "\n")
fmt(v::AbstractVector) = join([@sprintf("%10.5f", x) for x in v], " ")

"""
    save_results(rec, dir, tag)

Writes `tag_summary.txt` (final point) and `tag_history.csv` (certificate and
its components, prices and excess supplies by iteration) in `dir`.
"""
function save_results(rec, dir, tag)
    mkpath(dir)
    p = rec[:prices]; L, S = size(p)
    open(joinpath(dir, tag * "_summary.txt"), "w") do io
        println(io, "iterations (augmented Walrasian): ", rec[:iterations])
        @printf(io, "certificate after augmented Walrasian: %.3e\n", rec[:awlrs_certificate])
        @printf(io, "certificate (final): %.3e\n", rec[:certificate])
        println(io, "agent problems solved: ", rec[:nsolves])
        println(io, "portfolio bound active: ", rec[:bound_active])
        println(io, "rank of M1 at the final prices: ", rec[:rank_M1], " (contracts: ", length(rec[:q]), ")")
        println(io, "modified prices p~ (rows goods, columns nodes):\n", fmt(p))
        sig = p[1, :]
        println(io, "prices p = p~[:,s] / p~[1,s]:\n", fmt(p ./ sig'))
        println(io, "state prices sigma: ", fmt(sig[2:end]))
        println(io, "contract prices q: ", fmt(rec[:q]))
        @printf(io, "interest rate 1/q_1 - 1: %.4f\n", 1 / rec[:q][1] - 1)
        haskey(rec, :bliss) && @printf(io, "largest utility index at the solution: %.4f (bliss level K = %.4g)\n", rec[:bliss]...)
        println(io, "ES:\n", fmt(rec[:ES]))
        println(io, "Z: ", fmt(rec[:Z]))
        for (i, s) in enumerate(rec[:allocations])
            println(io, "agent $i: net portfolio z = ", fmt(s.z))
            println(io, "agent $i: consumption c =\n", fmt(s.c))
            any(!iszero, s.w) && println(io, "agent $i: retention w =\n", fmt(s.w))
        end
    end
    open(joinpath(dir, tag * "_history.csv"), "w") do io
        println(io, join(vcat(["it", "certificate", "ES_neg", "complementarity", "Z_inf", "r"],
                              ["p_$(l)_$(s)" for s in 1:S for l in 1:L],
                              ["ES_$(l)_$(s)" for s in 1:S for l in 1:L],
                              ["Z_$(j)" for j in eachindex(rec[:Z])]), ","))
        for h in rec[:history]
            println(io, join(vcat(h.it, h.cert, h.neg, h.comp, h.clr, h.r, vec(h.p), vec(h.ES), h.Z), ","))
        end
    end
    return dir
end

end # module
