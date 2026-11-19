# Changes (Julia port, 11 Nov 2026)

## Relative to the corrected Python code (September 2026)
- Language: Julia with JuMP and Ipopt (in process), NLopt (BOBYQA). Each agent
  problem is built once; prices are JuMP parameters; solves are warm started.
- Dual domain: fixed box G = [0, gmax]^(L S) including the numeraire coordinate
  (bounded dual domain of the revised Section 3). In Python the bound Mp grew
  along the iterations and was also used as the price bound; the price box is
  now separate (`pmin`, `pmax0`, `pmax1`).
- ES follows the paper: contract deliveries D_s z are included in the
  second-stage excess supply. The Python code left them out.
- Contract penalty weight rho: 1 to 3 (was 1e-2 to 1e-1). With deliveries inside
  ES, a small rho lets a contract imbalance pay for goods imbalances.
- Phase II tolerances: ftol 1e-6, xtol 1e-6, 40 evaluations per price unknown
  (was ftol 1e-2, xtol 1e-4). With the loose values BOBYQA stopped after its
  initial model and the iterations stalled.
- Minimum-norm portfolio selection is applied whenever M1(p) is rank deficient
  (was: only with more contracts than scenarios).
- Home production variables removed (not in the paper).
- Refinement: own Levenberg-Marquardt with forward differences (was SciPy
  least_squares).
- Parallel agent solves dropped (the Python parallel path ignored rho_w).
- Output: text summary and CSV history. `make_figures.jl` draws the figures
  from the CSV files with pgfplots.

## Example 4.3 (changed to fit the hypotheses of the revised paper)
- Utility: sum_s pi_s [ (prod c^b)^(1/2) + rho_w (prod w^b)^(1/2) ], rho_w = 0.01.
  Concave and nondecreasing everywhere. The quadratic bliss utility with
  K = 5.7 is not nondecreasing on the choice box (aggregate supplies reach 14),
  and with K raised above the box bound the iterations did not converge.
- Contracts: 2 real contracts with 3 scenarios (incomplete, M1 of full column
  rank): the bond, and a contract delivering good 2 in amounts (1.5, 0.5, 1.0).
  A plain forward on good 2 made M1 nearly rank deficient (portfolios of size
  10 to 60) and the refinement did not converge.
- Retention: the numeraire is carried to the second stage (A_s = e1 e1'), and
  retained goods enter the utility. With A_s = I (all goods carried over) the
  iterations did not converge.
- No transaction costs (D_0 = 0). Endowments: base vectors of the paper, scaled
  by 1 + 0.2 u in scenario 1 and 1 - 0.2 u in scenario 3, u uniform, seed 10.
- Start: goods priced like the numeraire, state prices equal to the
  probabilities, 10% random perturbation.

## Examples 4.1, 4.2
- Data unchanged. Example 4.2.1: aggregate endowment of good 2 reaches 7 > K, so
  the utility is not nondecreasing on the whole box; the largest utility index
  at the solution is 3.10 < K (reported in the summary).
- Example 4.2.2 does not converge (see the entry of 19 Nov 2026).

## 16 Nov 2026
- `solve_agent!`: on a solver failure, one retry from a cold start (endowment, zero portfolio). It is triggered only on failure, so runs without failures are unchanged.
- `experiments/lm_vs_awlrs.jl`: comparison of Levenberg-Marquardt started at the augmented Walrasian output and at the initial prices (10 starts per example).

## 19 Nov 2026 (release for the revised paper)
- The retry of 16 Nov is now an option, `retry` in `Economy` (default false).
  `run_all.jl` runs without it and reproduces `results/` of 11 Nov exactly
  (checked for the four examples); `experiments/lm_vs_awlrs.jl` runs with it,
  as on 16 Nov. With the retry, Example 4.3 in `run_all.jl` follows a different
  path (a failed agent solve at a Phase II trial price is rescued instead of
  penalized) and ends at another point of the same set of equilibrium prices,
  with the same (p, q).
- Comments corrected: Example 4.2.2 was described as having equilibrium prices
  at which M1(p) drops rank. At an equilibrium the relative prices of good 2 in
  the three scenarios cannot coincide (Section 4.2.2 of the paper); the computed
  prices are close to each other and M1(p) is close to losing rank.
- Comments in `examples/ex4.jl` updated (retention of the numeraire only,
  power utility); reference to equation (7) of the paper in `contract_prices`;
  correspondence between `r` and r^nu of the paper in `Options`.
- `scripts/checks.jl`: checks of Section 4 and Table 4; output in `results/checks.txt`.
- `experiments/reproduce.sh`; `experiments/summarize.py` reads its results
  relative to its own folder. The Example 4.3 files `*_seeds1-10.csv` were
  renamed `*_seeds1-7.csv`, the seeds they contain.
- README rewritten for the revised paper.
