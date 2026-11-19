# Figures of Section 4 from results/*_history.csv, typeset by pgfplots so the
# fonts are those of the paper (Computer Modern). Needs pdflatex with pgfplots.
#
#   julia make_figures.jl            # black-and-white and colour versions
#
# Output in figures/ (one PDF per figure and version, suffix _bw or _color):
#   convergence_all    certificate and its three components, four examples
#   <tag>_convergence  the same for one example
#   <tag>_prices       modified prices by node (small multiples)
#   <tag>_excess       excess supply by node (small multiples), iterations >= 1
#
# Both versions use the same line patterns and marks; the colour version only
# adds colour, so it stays readable when printed in black and white.

const ROOT = @__DIR__
const RES = joinpath(ROOT, "results")
const FIG = joinpath(ROOT, "figures")
const EXAMPLES = [("ex1_bde", "Example 4.1"), ("ex2_sch2", "Example 4.2.1"),
                  ("ex2_sch3", "Example 4.2.2"), ("ex4", "Example 4.3")]

# line patterns for the goods (index 1 = numeraire) and for the certificate parts
const PATTERN = ["solid", "densely dashed", "densely dotted", "dashdotted",
                 "dashdotdotted", "loosely dashed"]
const MARK = ["*", "square*", "triangle*", "o", "square", "triangle"]
const GRAY = ["black", "black!75", "black!55", "black", "black!75", "black!55"]
const COLOR = ["black", "cbblue", "cborange", "cbgreen", "cbred", "cbpurple"]

preamble() = raw"""
\documentclass[border=2pt]{standalone}
\usepackage{amsmath}
\usepackage{pgfplots}
\usepgfplotslibrary{groupplots}
\pgfplotsset{compat=1.16}
\definecolor{cbblue}{RGB}{0,114,178}
\definecolor{cborange}{RGB}{213,94,0}
\definecolor{cbgreen}{RGB}{0,158,115}
\definecolor{cbred}{RGB}{204,121,167}
\definecolor{cbpurple}{RGB}{86,180,233}
\pgfplotsset{
  every axis/.append style={font=\footnotesize, tick label style={font=\scriptsize},
    label style={font=\footnotesize}, title style={font=\footnotesize, yshift=-3pt},
    legend style={font=\scriptsize, draw=none, fill=none, legend cell align=left},
    axis line style={line width=0.4pt}, tick style={line width=0.4pt, black},
    grid=major, grid style={line width=0.2pt, black!15}},
  every axis plot/.append style={line width=0.7pt, mark size=1.1pt},
}
\begin{document}
"""

function shape(tag)
    cols = split(readline(joinpath(RES, tag * "_history.csv")), ",")
    pc = [parse.(Int, split(c, "_")[2:3]) for c in cols if startswith(c, "p_")]
    nrows = countlines(joinpath(RES, tag * "_history.csv")) - 1
    return maximum(first, pc), maximum(last, pc), nrows - 1   # L, S, last "it"
end

style(k, bw) = string(PATTERN[k], ", ", bw ? GRAY[k] : COLOR[k])

"Axis content: certificate and components on a log scale; the last point is the refinement."
function convergence_axis(tag, bw; legend = false)
    L, S, last = shape(tag)
    file = "../results/$(tag)_history.csv"
    series = [("certificate", "certificate", 1, "line width=1.1pt"),
              ("ES_neg", raw"$\|ES^-\|_\infty$", 2, ""),
              ("complementarity", raw"$\max_{\xi,l}\tilde p_{\xi,l}|ES_{\xi,l}|$", 3, ""),
              ("Z_inf", raw"$\|Z\|_\infty$", 4, "")]
    io = IOBuffer()
    for (col, lab, k, extra) in series
        # augmented Walrasian iterations (all rows but the last)
        println(io, "\\addplot[$(style(k, bw)), $extra, unbounded coords=jump] ",
                "table[col sep=comma, x=it, y=$col, restrict x to domain=0:$(last - 1)] {$file};")
        legend && println(io, "\\addlegendentry{$lab}")
        # refinement step (last row), joined by a thin line
        println(io, "\\addplot[$(bw ? GRAY[k] : COLOR[k]), line width=0.3pt, mark=$(MARK[k]), ",
                "mark size=1.4pt, forget plot, unbounded coords=jump] ",
                "table[col sep=comma, x=it, y=$col, restrict x to domain=$(last - 1):$last] {$file};")
    end
    println(io, "\\draw[black!50, line width=0.3pt, densely dotted] (axis cs:$(last - 0.5),1e-9) -- (axis cs:$(last - 0.5),1e3);")
    return String(take!(io)), last
end

function fig_convergence(tag, title, bw)
    body, last = convergence_axis(tag, bw; legend = true)
    """
    \\begin{tikzpicture}
    \\begin{semilogyaxis}[width=8.4cm, height=5.6cm, xlabel={iteration \$\\nu\$}, title={$title},
      xmin=0, xmax=$(last + 0.5), ymin=1e-8, ymax=1e3, clip=true,
      ytick={1e-8,1e-6,1e-4,1e-2,1e0,1e2}, legend pos=south west,
      extra x ticks={$last}, extra x tick labels={ref.}, extra x tick style={grid=none, tick label style={yshift=-7pt}}]
    $body
    \\end{semilogyaxis}
    \\end{tikzpicture}
    """
end

function fig_convergence_all(bw)
    io = IOBuffer()
    println(io, raw"\begin{tikzpicture}")
    println(io, raw"\begin{groupplot}[group style={group size=2 by 2, horizontal sep=0.9cm, vertical sep=1.25cm,",
            raw" ylabels at=edge left, yticklabels at=edge left}, ymode=log,",
            raw" width=6.6cm, height=4.6cm, ymin=1e-8, ymax=1e3, ytick={1e-8,1e-6,1e-4,1e-2,1e0,1e2}, xmin=0]")
    for (n, (tag, title)) in enumerate(EXAMPLES)
        body, last = convergence_axis(tag, bw; legend = n == 1)
        opts = "title={$title}, xmax=$(last + 0.5), extra x ticks={$last}, extra x tick labels={ref.}, " *
               "extra x tick style={grid=none, tick label style={yshift=-7pt}}"
        n >= 3 && (opts *= ", xlabel={iteration \$\\nu\$}")
        n == 1 && (opts *= ", legend to name=convlegend, legend columns=4, legend style={/tikz/every even column/.append style={column sep=6pt}}")
        println(io, "\\nextgroupplot[$opts]")
        print(io, body)
    end
    println(io, raw"\end{groupplot}")
    println(io, raw"\node[anchor=north] at ($(group c1r2.south)!0.5!(group c2r2.south)+(0,-0.95cm)$) {\pgfplotslegendfromname{convlegend}};")
    println(io, raw"\end{tikzpicture}")
    String(take!(io))
end

"Small multiples by node of the columns prefix_l_s; `from` is the first iteration shown."
function fig_panels(tag, prefix, ylabel, bw; from = 0, zero_line = false)
    L, S, last = shape(tag)
    file = "../results/$(tag)_history.csv"
    io = IOBuffer()
    println(io, raw"\begin{tikzpicture}")
    println(io, "\\begin{groupplot}[group style={group size=$S by 1, horizontal sep=0.75cm}, ",
            "width=$(round(13.6 / S + 1.3, digits = 2))cm, height=4.6cm, xmin=$from, xmax=$last, xlabel={iteration \$\\nu\$}]")
    for s in 1:S
        opts = "title={" * (s == 1 ? "first stage" : "scenario \$\\xi=$(s - 1)\$") * "}"
        s == 1 && (opts *= ", ylabel={$ylabel}, legend to name=goodlegend, legend columns=$L, " *
                           "legend style={/tikz/every even column/.append style={column sep=6pt}}")
        println(io, "\\nextgroupplot[$opts]")
        zero_line && println(io, "\\addplot[black!40, line width=0.3pt, forget plot, domain=$from:$last, samples=2] {0};")
        for l in 1:L
            println(io, "\\addplot[$(style(l, bw)), mark=$(MARK[l]), mark repeat=$(max(1, last ÷ 6)), mark phase=$(l)] ",
                    "table[col sep=comma, x=it, y=$(prefix)_$(l)_$(s), restrict x to domain=$from:$last] {$file};")
            s == 1 && println(io, "\\addlegendentry{good \$l=$(l - 1)\$}")
        end
    end
    println(io, raw"\end{groupplot}")
    println(io, "\\node[anchor=north] at (\$(group c1r1.south west)!0.5!(group c$(S)r1.south east)+(0,-0.85cm)\$) {\\pgfplotslegendfromname{goodlegend}};")
    println(io, raw"\end{tikzpicture}")
    String(take!(io))
end

function build(name, body)
    mkpath(FIG)
    tex = joinpath(FIG, name * ".tex")
    write(tex, preamble() * "\\usetikzlibrary{calc}\n" * body * "\n\\end{document}\n")
    for _ in 1:2   # second pass resolves the legend references
        run(pipeline(Cmd(`pdflatex -interaction=nonstopmode -halt-on-error $(name).tex`; dir = FIG),
                     stdout = devnull))
    end
    foreach(e -> rm(joinpath(FIG, name * e); force = true), (".aux", ".log"))
    println("figures/", name, ".pdf")
end

for (bw, suffix) in ((true, "_bw"), (false, "_color"))
    build("convergence_all" * suffix, fig_convergence_all(bw))
    for (tag, title) in EXAMPLES
        isfile(joinpath(RES, tag * "_history.csv")) || continue
        build(tag * "_convergence" * suffix, fig_convergence(tag, title, bw))
        build(tag * "_prices" * suffix, fig_panels(tag, "p", raw"modified prices $\tilde p_{\xi,l}$", bw))
        build(tag * "_excess" * suffix, fig_panels(tag, "ES", raw"excess supply $ES_{\xi,l}$", bw; from = 1, zero_line = true))
    end
end
