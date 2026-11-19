import csv, glob, os, statistics as st, math
rows = []
for f in sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', '*.csv'))):
    rows += [r for r in csv.DictReader(open(f)) if r['status'] != 'start not evaluable']
order = ['ex1_bde', 'ex2_sch2', 'ex2_sch3', 'ex4']
name = {'ex1_bde': '4.1', 'ex2_sch2': '4.2.1', 'ex2_sch3': '4.2.2', 'ex4': '4.3'}
def fl(x):
    try: return float(x)
    except: return math.nan
def ok(r):  # equilibrium of the original economy: kappa <= 1e-6 and no contract bound active
    return r['status'] == 'ok' and fl(r['kappa_final']) <= 1e-6 and r['bound_active'] == 'false'
def okb(r): # equilibrium of the bounded economy only
    return r['status'] == 'ok' and fl(r['kappa_final']) <= 1e-6 and r['bound_active'] == 'true'
methods = ['AW+LM', 'AW10+LM', 'AW5+LM', 'LM']
out = []
out.append("| Ex. | method | starts | success | bound active | median final kappa | max final kappa | median agent solves | min | max | share in AW phase | median time (s) |")
out.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
for tag in order:
    for m in methods:
        R = sorted([r for r in rows if r['tag'] == tag and r['method'] == m], key=lambda r: int(r['seed']))
        if not R: continue
        k = [fl(r['kappa_final']) for r in R]; kf = [x for x in k if not math.isnan(x)]
        s = [int(r['solves_total']) for r in R]; aw = [int(r['solves_aw']) for r in R]
        out.append(f"| {name[tag]} | {m} | {len(R)} | {sum(map(ok,R))}/{len(R)} | {sum(map(okb,R))} | {st.median(kf):.1e} | {max(kf):.1e} | {int(st.median(s))} | {min(s)} | {max(s)} | {sum(aw)/sum(s):.2f} | {st.median(fl(r['time_s']) for r in R):.1f} |")
print('\n'.join(out))
print()
# LM with restarts: from start s, run LM on starts s, s+1, ... (cyclically) until success
print("LM with restarts (start s, then s+1, ... cyclically, until success):")
for tag in order:
    R = {int(r['seed']): r for r in rows if r['tag'] == tag and r['method'] == 'LM'}
    if not any(ok(r) for r in R.values()):
        print(f"  {name[tag]}: no LM success"); continue
    seeds = sorted(R); costs = []; tries = []
    for i in range(len(seeds)):
        c = 0
        for j in range(len(seeds)):
            r = R[seeds[(i + j) % len(seeds)]]; c += int(r['solves_total'])
            if ok(r): break
        costs.append(c); tries.append(j + 1)
    print(f"  {name[tag]}: success {len(seeds)}/{len(seeds)}, median solves {int(st.median(costs))}, max {max(costs)}, max starts used {max(tries)}")
print()
print("Final (q0, q1) of the successful runs, by example:")
for tag in order:
    qs = sorted({(r['q0'], r['q1']) for r in rows if r['tag'] == tag and ok(r)})
    print(f"  {name[tag]}: {qs}")
for tag in order:
    awk = [fl(r['kappa_aw']) for r in rows if r['tag'] == tag and r['method'] == 'AW+LM']
    awk = [x for x in awk if not math.isnan(x)]
    if awk: print(f"  {name[tag]}: kappa after the AW iterations: median {st.median(awk):.1e}, range [{min(awk):.1e}, {max(awk):.1e}]")
