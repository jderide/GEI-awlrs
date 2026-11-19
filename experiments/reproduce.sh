#!/bin/sh
# Commands that produced experiments/results/ (Section 4.4, Table 7), from the
# repository root. About 2 hours in total on the machine used (plus compilation); Example 4.3 takes most of it.
# The Example 4.3 runs were done in two parts: seeds 1 to 7, then 8 to 11.
# At seed 8 an agent's problem cannot be solved at the starting prices, so the
# start is excluded and seed 11 replaces it (Section 4.4).
set -e
for tag in ex1_bde ex2_sch2 ex2_sch3; do
  julia --project=. experiments/lm_vs_awlrs.jl $tag 1 10        # AW+LM and LM
  julia --project=. experiments/lm_vs_awlrs.jl $tag 1 10 5      # AW5+LM
  julia --project=. experiments/lm_vs_awlrs.jl $tag 1 10 10     # AW10+LM
done
julia --project=. experiments/lm_vs_awlrs.jl ex4 1 7
julia --project=. experiments/lm_vs_awlrs.jl ex4 8 11
julia --project=. experiments/lm_vs_awlrs.jl ex4 1 7 5
julia --project=. experiments/lm_vs_awlrs.jl ex4 8 11 5
python3 experiments/summarize.py > experiments/results/summary.txt
