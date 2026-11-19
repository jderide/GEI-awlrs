# GEI-awlrs: augmented Walrasian for GEI economies with retention

Code, data and outputs for

> J. Deride, *Solving equilibrium problems in economies with incomplete financial markets and retention*, Set-Valued and Variational Analysis (revised version, November 2026).

The code computes equilibrium prices of the modified economy of the paper by the augmented Walrasian iterations of Section 3, followed by a Levenberg–Marquardt refinement, and reproduces Section 4.

## Contents

| Path | Content |
|---|---|
| `src/GEIawlrs.jl` | model, agents' problems (JuMP, Ipopt), excess supply, augmented Walrasian iterations (Phase II by BOBYQA, NLopt), refinement, output |
| `examples/ex1_bde.jl` | Example 4.1 (Brown, DeMarzo, Eaves) |
| `examples/ex2_sch2.jl` | Example 4.2.1 (Schmedders, endowment variation) |
| `examples/ex2_sch3.jl` | Example 4.2.2 (Schmedders, one agent of each type) |
| `examples/ex4.jl` | Example 4.3 (5 agents, 6 goods, retention, 2 real contracts) |
| `run_all.jl` | runs the four examples with seed 10; writes `results/` |
| `results/` | `table.txt` (Table 3), per example `*_summary.txt` (final point) and `*_history.csv` (iterations), and `checks.txt` (output of `scripts/checks.jl`) |
| `make_figures.jl` | Figures 2 to 5 from `results/*_history.csv` (pgfplots, needs pdflatex); black-and-white and colour versions in `figures/` |
| `scripts/checks.jl` | statements of Section 4 evaluated at the final prices (ranks, bounds, positions, Example 4.3 endowment factors) and Table 4 |
| `experiments/lm_vs_awlrs.jl` | Section 4.4: refinement started after the iterations (AW+LM) and from the starting prices (LM) |
| `experiments/reproduce.sh` | the commands that produced `experiments/results/` |
| `experiments/summarize.py` | Table 7 from `experiments/results/*.csv`; output in `experiments/results/summary.txt` |
| `CHANGES.md` | history of the code |

## Install and run

The stored outputs in `results/` and `experiments/results/` were computed on Linux x86_64 with Julia 1.13.3; the package versions (JuMP, Ipopt, NLopt) are fixed in `Manifest.toml`.
On another platform or Julia version the iterates differ from the stored ones by rounding, and the differences grow along the iterations: $\kappa$ along the iterations, the number of agent solves and the modified prices $\tilde p$ change.
The prices $(p,q)$ do not: a run on macOS gave the same $(p,q)$ to five digits in Examples 4.1, 4.2.1 and 4.3.

    julia --project=. -e 'using Pkg; Pkg.instantiate()'
    julia --project=. run_all.jl              # Table 3, results/ (about 8 minutes)
    julia --project=. scripts/checks.jl > results/checks.txt   # after run_all.jl
    julia make_figures.jl                     # figures/, after run_all.jl
    sh experiments/reproduce.sh               # Table 7 (about 2 hours)
    julia --project=. examples/ex1_bde.jl     # one example, with the iteration log

## Where each item of Section 4 comes from

| Paper | Source |
|---|---|
| Table 2 (parameters) | `options_bde()` in `examples/ex1_bde.jl` (Examples 4.1, 4.2.1, 4.2.2), `options_ex4()` in `examples/ex4.jl`, defaults in `Options` (`src/GEIawlrs.jl`); see the correspondence below |
| Table 3, prices and interest rates of Examples 4.1, 4.2.1, 4.3 | `run_all.jl`, `results/table.txt`, `results/*_summary.txt` |
| Figures 2 to 5 | `make_figures.jl` |
| rank of $M_1$, bounds $m$ and $M$, contract positions, index 3.10 in Example 4.2.1, relative prices and positions in Example 4.2.2, factors $a$, $b$ in Example 4.3 | `scripts/checks.jl` |
| Table 4 | `scripts/checks.jl` |
| Table 7 and the text of Section 4.4 | `experiments/`, `experiments/results/summary.txt` |

## Notation

| Paper | Code |
|---|---|
| goods $l=0,\dots,L$, numéraire $l=0$ | rows `1:L` of the price matrix, numéraire row 1 |
| stage and states $\xi=0,\dots,\Xi$ | columns `1:S` (node 1 is the first stage) |
| contracts $j=0,\dots,J$, bond $j=0$ | contracts `1:J`, bond contract 1 |
| modified prices $\tilde p$ | `p` (with `p[1,1] = 1`) |
| $q$ by (7) | `contract_prices` |
| $M_1(\tilde p)$ | `payoff_values` |
| $\mathrm{ES}$, $Z$ | `excess` |
| $\kappa$ | `certificate` (the code calls $\kappa$ the certificate) |
| bound $m$, multiple 1.03 | `ub` in `Economy` |
| bound $M = 100$ | `zbound` |
| box $B$ | `pmin`, `pmax0`, `pmax1` |
| $\bar g^\nu$ | `gmax` |
| $\rho$ | `rho` |
| $r^\nu$ | `1/r`: the code uses $r$ in place of $1/r^\nu$; `r0 = 0.25`, `r1 = 5` give $r^\nu$ from 4 to 0.2, `r1 = 10` gives 0.1 |
| step $\delta$ | `damping` |
| BOBYQA evaluations per price unknown | `maxeval` |
| tolerance $\varepsilon$, iteration limit | `tol`, `maxit` |

The agents' problems are solved with Ipopt (`tol = 1e-10`, `acceptable_tol = 1e-8`, `max_iter = 500`).
In `run_all.jl` (Table 3) a failure on an agent's problem is penalized: Phase II assigns a very low value and the refinement a large residual to that price.
In `experiments/` (Section 4.4) the problem is first solved once more from the endowment and a zero portfolio, and both solves are counted (`retry = true` in `Economy`).
Note on the revised manuscript of 19 Nov 2026: Section 4 describes this retry for all runs; the runs of Table 3 and Figures 2 to 5 were made without it, as above.
With the retry, Example 4.3 in `run_all.jl` ends at another point of the same set of equilibrium prices, with the same $(p,q)$, after a different number of iterations and agent solves.
Phase II uses BOBYQA with `ftol_rel = xtol_rel = 1e-6`.
When $M_1(\tilde p)$ does not have full column rank, each agent's net position is replaced by the net position of minimum norm with the same payoffs.

## Results (seed 10, `results/table.txt`)

| Example | iterations | $\kappa$, last iterate | $\kappa$, refined | interest rate |
|---|---|---|---|---|
| 4.1 | 50 | 3.4e-2 | 6.6e-8 | 8.65% |
| 4.2.1 | 50 | 2.5e-2 | 8.0e-8 | 1.73% |
| 4.2.2 | 50 | 4.8e-1 | 1.8e-1 | no equilibrium found |
| 4.3 | 30 | 2.4e-1 | 3.4e-8 | 2.50% |

In Example 4.2.2 neither the iterations nor the refinement reach an equilibrium: the computed relative prices of good 2 in the three scenarios are close to each other, $M_1(\tilde p)$ is close to losing rank, and the contract positions are large (Section 4.2.2 of the paper).

## License

GNU General Public License v3.0 (`LICENSE`).
