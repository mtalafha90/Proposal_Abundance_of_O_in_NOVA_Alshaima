# Results

Every number below comes from the runs in `results/`, computed with NucNetPy.
Figures are in `figures/`; the LaTeX versions of the tables are in
`results/tables/`.

The trajectory model uses
`data/trajectories/iliadis2002_S1_synthetic_benchmark.txt`, which starts at
$T_9 = 0.070$, $\rho = 2.200\times10^4$ g cm$^{-3}$, peaks at $T_9 = 0.418$ at
$t = 100.0$ s where $\rho = 4.000\times10^3$ g cm$^{-3}$, and is followed to the
end of the file at $t = 3000$ s.

> **What this profile is.** It is a *literature-constrained synthetic
> benchmark*, **not** a hydrodynamic trajectory. Only the peak temperature
> (0.418 GK) and the model identity — ONe white dwarf, 1.35 M☉, 50% mixing,
> model S1 of Iliadis et al. (2002, ApJS 142, 105) — are taken from the
> literature. Iliadis et al. do not publish the time series, so the *shape* is
> an analytic construction: a logistic rise to the peak at 100 s, density
> falling log-linearly with the same rise function, then stretched-exponential
> cooling ($\tau_T = 25$ s, $\alpha_T = 1.25$) and expansion ($\tau_\rho = 45$ s,
> $\alpha_\rho = 0.9$). The generator, its parameter summary and the upstream
> README are in `data/trajectories/provenance/`.
>
> **Consequence:** both thermodynamic histories in this study are now analytic.
> The comparison is between a simple exponential cooling and expansion
> prescription and a richer parameterisation, which is well posed, but nothing
> here is evidence of what a
> multi-zone hydrodynamic calculation would give. Do not cite these results as
> results of Iliadis et al. (2002) beyond the peak temperature anchor.


## 1. The two reference calculations

| Model | Thermodynamic behaviour | $R^{\rm initial}$ | $R^{\rm final}$ | $f_{\rm enh}$ | $t_{\rm fo}$ | $T_{9,\max}$ |
|---|---|---|---|---|---|---|
| Exponential | immediate cooling | 3.68e-3 | **0.460** | 125 | 0.23 s | 0.200 |
| Trajectory | delayed temperature peak | 3.68e-3 | **1.958** | 533 | 144.5 s | 0.418 |

Both agree with the values quoted in the proposal, which gives
$R^{\rm final} \approx 0.4$–$0.5$ for the exponential model and "of order a
few" for the trajectory model. The shapes of the curves also match the
proposal's figures: an early bump, a sharp dip, a long climb, an overshoot, and
a flat freeze-out plateau (`figures/fig07_ratio_comparison.png`).

The ratio peaks at 2.49 during the burning phase and settles at 1.958, so the
trajectory model enhances $^{15}$N-rich material by a factor of 533 over the
solar starting composition, against 125 for the exponential model.

**One point of bookkeeping.** The methodology chapter defines the ratio from
molar abundances,
$R = (Y_{15\rm N}+Y_{15\rm O})/(Y_{14\rm N}+Y_{14\rm O})$, which for the solar
starting composition gives $R^{\rm initial} = 3.68\times10^{-3}$. Chapter 3
quotes $3.9\times10^{-3}$, which is the same composition evaluated as a ratio of
*mass fractions*, $X(^{15}{\rm N})/X(^{14}{\rm N})$. The two differ by the
factor $14/15$. The molar definition is used throughout here, because that is
the one the methodology chapter gives; the enhancement factor is barely
affected either way.

## 2. Why the ratio moves when it does

The reaction flows (`fig08`, `fig09`) and the timescales (`fig10`, `fig11`)
explain the shape of the curves.

**The early dip is real, not numerical.** At the starting temperature the two
fastest reactions in the whole network are $^{18}{\rm O}({\rm p},\alpha)^{15}{\rm N}$
and $^{15}{\rm N}({\rm p},\alpha)^{12}{\rm C}$. The second destroys the solar
$^{15}$N far faster than anything replaces it, so $R$ falls by more than an
order of magnitude before the CNO cycle has begun to turn. This is also what
sets the initial nuclear energy generation, about $10^{15}$ erg g$^{-1}$
s$^{-1}$ (`fig14`, `fig15`).

**The climb is the hot-CNO cycle turning over.** $^{12}{\rm C}({\rm p},\gamma)$
feeds $^{13}$N, and from there the flow reaches $^{14}$O and $^{15}$O. Because
$^{14}$O and $^{15}$O are beta-unstable with half-lives of 71 s and 122 s, they
build up while the material is hot and then decay to $^{14}$N and $^{15}$N once
it is not. $^{15}$O is made faster than $^{14}$O, so the ratio ends far above
its starting value.

**Freeze-out: two criteria, and they are not the same one.** $t_{\rm fo}$ is
measured numerically — the first time after which $R$ stays within 1% of its
final value. The timescale crossing is a *separate* criterion that explains it.
They agree, but not as a single number, because the reactions that can still
move $R$ do not stop together:

| Reaction | crosses $\tau_T$ (trajectory) | (exponential) |
|---|---|---|
| $^{14}$N(p,γ) — feeds A=15 | 128.1 s (T₉ 0.138) | 0.14 s |
| $^{12}$C(p,γ) — feeds the cycle | 134.0 s (T₉ 0.104) | 0.23 s |
| $^{15}$N(p,α) — drains A=15 | 147.2 s (T₉ 0.055) | 0.59 s |
| **numerical $t_{\rm fo}$** | **144.5 s** | **0.23 s** |

In both models $t_{\rm fo}$ falls inside the interval the crossings span, close
to the last one — the reaction that keeps changing the ratio longest. That is
the honest relationship: the crossing interprets $t_{\rm fo}$, it does not
define it.

**How the thermodynamic timescales are computed.** $T_9$ and $\rho$ come from
the same thermo callable the solver uses, sampled at the output times; for the
trajectory that is linear interpolation between table rows. No smoothing,
anywhere. Derivatives are `numpy.gradient` — second-order central differences on
the non-uniform grid, second-order one-sided at the ends. Through the burning
episode the S1 table has 0.050 s spacing and the output grid 0.042 s, so the
output is marginally the finer of the two and consecutive output points can
fall inside one table interval, where the linear interpolant has constant
slope. That makes d$T_9$/d$t$ mildly stepped at the 0.05 s scale. It does not
affect the comparison used here, which is made on the cooling branch over tens
of seconds, but it is the reason $\tau_T$ should not be read at the resolution
of a single output point. (The S1 file also lists its peak row twice, at
t = 100 s; `nova_trajectory` drops the repeat, which leaves the interpolant
bit-identical and keeps any derivative of the table finite.) $\tau_T$ diverges where
d$T_9$/d$t$ = 0 (at the peak and at minor stationary points); it is not
regularised, is drawn with a ceiling in the figures, and the nuclear-versus-
thermodynamic comparison is only used on the cooling branch where d$T_9$/d$t$ is
large — which is where freeze-out happens anyway.

**The trajectory produces more $^{15}$N because it burns hotter, for longer.**
The peak temperature is 0.418 against 0.200, and the profile holds the material
above $T_9 = 0.2$ for 38.4 seconds and above $T_9 = 0.1$ for 64.2 seconds rather
than cooling away from the peak immediately. Those are two of the five
things that differ between the histories; section 6 separates them with matched
controls and finds the exposure time to be the largest single contribution.

## 3. What happens to oxygen

The proposal is about oxygen in novae, so it is worth recording the elemental
outcome as well as the isotopic ratio. Total elemental mass fractions, start to
finish:

(Initial values are the solar abundances after the renormalisation described in
section 8 — removing species outside the network's charge range raises them by
0.32% for `Z ≤ 10` — so they are what the runs actually started from.)

| Element | Initial | Exponential | Trajectory ($Z\leq10$) | Trajectory ($Z\leq30$) |
|---|---|---|---|---|
| C | 2.31e-3 | 1.11e-4 (×0.05) | 5.01e-4 (×0.22) | 4.98e-4 (×0.22) |
| N | 8.03e-4 | 2.27e-3 (×2.8) | 8.36e-3 (×10.4) | 8.32e-3 (×10.4) |
| O | 5.81e-3 | 6.78e-3 (×1.2) | 1.22e-5 (×2.1e-3) | 3.32e-7 (×5.7e-5) |

**The two thermodynamic histories do opposite things to oxygen.** The
exponential model, which never gets above $T_9 = 0.2$, leaves oxygen slightly
enhanced: it burns carbon into nitrogen and a little oxygen, and then freezes.
The trajectory model, which reaches $T_9 = 0.418$, destroys oxygen by a factor
of 475 and converts it into nitrogen, which ends up ten times its solar value.
$^{16}$O falls from a mass fraction of 5.77e-3 to 6.49e-6 (the elemental total
of 5.79e-3 includes $^{17}$O and $^{18}$O).

The path is the hot-CNO one:
$^{16}{\rm O}({\rm p},\gamma)^{17}{\rm F}(\beta^+)^{17}{\rm O}({\rm p},\alpha)^{14}{\rm N}$.
Once the temperature is high enough for $^{16}$O to capture a proton faster
than the envelope cools, oxygen is a way station on the road to nitrogen
rather than an endpoint. This is the well-known signature of nova ejecta:
nitrogen-rich, carbon- and oxygen-poor relative to solar.

So the hot, prolonged history that raises $R_{15/14}$ by a factor of seven over
the exponential case also decides whether oxygen survives at all. Any
measurement of oxygen in nova ejecta is therefore a strong constraint on the
peak thermodynamic conditions — the effect on it is three orders of magnitude
against a factor of seven for the isotopic ratio. Note the caveat in section 7:
elemental oxygen is less network-converged than $R_{15/14}$.

## 4. Steady flow

Testing for steady flow needs **every link of the closed cycle**, not a
selection of proton captures: in the beta-limited regime the two decays are the
slow steps, so a test built only from the captures cannot show whether the cycle
circulates uniformly. All six links of

    12C(p,g)13N(p,g)14O(b+)14N(p,g)15O(b+)15N(p,a)12C

are compared, summarised by `D = max_j |log10 Q_j|` — the largest deviation of
any link from the reference flow. Taking the maximum, not an average, is what
makes it a whole-cycle test.

The window matters too. Ratios are meaningless where the denominator is
collapsing, so the analysis is restricted to where the **cycle throughput** (the
smallest flow around the loop) exceeds a tenth of its maximum: t = 59.0–126.0 s,
T₉ = 0.098–0.447. Using the reference flow alone would be worse — 14N(p,γ)
spikes sharply just before the peak while its target is being consumed, which is
the 14N reservoir emptying, not circulation.

**Result: the six principal links approach a common flow for a short interval
after the temperature maximum.**

| Criterion | Value |
|---|---|
| closest approach | D = 0.013 dex — all six links within **3%** |
| when | t = 104.12 s, 0.23 s after the T₉ maximum |
| all six within factor 1.2 | 0.43 s |
| all six within factor 2 | 2.73 s |
| all six within factor 3 | 13.6 s |

Resolved at 0.04 s output spacing across the burning episode. Before the peak
D ≈ 0.7–0.8 dex, rising to 2.2 dex during the 14N(p,γ) transient; afterwards it
grows back to 1.3 dex. So the cycle enters this condition abruptly, as the proton
captures overtake the decays, and leaves it gradually as the material cools.

### 4a. Side branches: the agreement is not steady flow

Equal flows around a loop are necessary for steady flow but not sufficient — a
nuclide can also be fed or drained by reactions outside the loop. Each
intermediate is therefore tested separately against the **whole** 68-nuclide
network, using

    f_side,i = sum|F_side,i| / (sum|F_principal,i| + sum|F_side,i|)

evaluated from the stored abundances of all 68 species (`traj_ref_abundances.csv`),
over the 2.73 s interval t = 103.86–106.59 s in which D < log10(2).

| nuclide | median f_side | max f_side |
|---|---|---|
| ¹²C | 7.2e-04 | 3.9e-02 |
| ¹³N | 1.5e-05 | 2.1e-05 |
| ¹⁴O | 1.1e-03 | 2.6e-02 |
| ¹⁴N | 4.3e-05 | 2.1e-01 |
| **¹⁵O** | **9.6e-01** | **9.8e-01** |
| ¹⁵N | 5.0e-05 | 5.2e-05 |

For five of the six, side flows are generally small over the interval, though
¹⁴N reaches 0.21 near its edge; their medians lie between 1.5e-05 and 1.1e-03,
so for most of the interval they are fed and drained by the cycle alone.
**¹⁵O is different.**
At t = 104.12 s (T₉ = 0.3954, ρ = 1987 g/cc):

| reaction | flow | share of ¹⁵O turnover |
|---|---|---|
| ¹⁸F(p,α)¹⁵O | 1.222e-04 | 98.1% |
| ¹⁵O(β⁺)¹⁵N | 1.151e-06 | 0.92% |
| ¹⁴N(p,γ)¹⁵O | 1.132e-06 | 0.91% |

¹⁸F(p,α)¹⁵O runs at **108×** the reference loop flow. Its source is
¹⁸Ne(β⁺)¹⁸F at 1.223e-04, which supplies >99.9% of the ¹⁸F; the alternative
¹⁷O(p,γ)¹⁸F contributes 9.7e-11, six orders of magnitude less. The path is the
CNO-II/III chain

    16O(p,g)17F(p,g)18Ne(b+)18F(p,a)15O

and along this trajectory it is the dominant supply of ¹⁵O, not a minor branch.

Consequently ¹⁵O is **not in balance**: production exceeds destruction by
1.21e-04, and X(¹⁵O) duly rises from 2.63e-03 at the temperature maximum to
3.04e-03 at t = 104.12 s, a rate that matches the imbalance to within a few per
cent. The closed cycle is therefore not closed at ¹⁵O, and the agreement of the
six links is **not** by itself evidence of a circulating steady state.

The three side-path flows are now written by the solver itself
(`F_f17_pg_ne18`, `F_ne18_bd_f18`, `F_f18_pa_o15`, `F_o17_pg_f18` in the flows
CSV), so none of the above rests on post-processing. ¹⁸F(p,α) exceeds
14N(p,γ) from just after the temperature maximum until t = 115.6 s, and
¹⁸Ne(β⁺) tracks ¹⁸F(p,α) to within 1% across the hot phase — ¹⁸F is destroyed
as fast as it is made.

**This is a property of the trajectory, not of the network.** In the
exponential model (T₉ peak 0.200) the dispersion never falls below 1.35 dex, so
the factor-2 condition is never met at all; at its closest approach
f_side(¹⁵O) = 1.2e-03, three orders of magnitude smaller. Building the ¹⁸Ne
reservoir via ¹⁶O(p,γ)¹⁷F(p,γ)¹⁸Ne needs the higher temperature the trajectory
reaches.

**The limiting step changes with temperature**, so the classical CNO
bottleneck does not hold throughout. The proton-capture timescale of
$^{14}$N is 244 s at the starting conditions, longer than the beta lifetimes
of $^{14}$O (101.9 s) and $^{15}$O (176.0 s), so $^{14}$N(p,$\gamma$) limits the
flow during the quiescent phase. It drops below the $^{15}$O lifetime at
$T_9 \simeq 0.10$ (t $\simeq$ 74 s) and reaches $3.9\times10^{-5}$ s at the
temperature maximum — six orders of magnitude faster. Through the whole burning
episode the cycle is therefore **beta-limited**.

**There are three waiting points, and ¹⁵O is not the largest.** At the
temperature maximum (t = 103.907 s, T₉ = 0.4467), 96% of the total abundance in
the hot-CNO species pool sits on three beta-unstable nuclei. The pool is all 40
C, N, O and F isotopes of the network (⁹⁻¹⁶C, ¹¹⁻¹⁸N, ¹²⁻²¹O, ¹⁴⁻²³F) plus the
proton-rich ¹⁶⁻¹⁹Ne; ²⁰⁻²⁵Ne are excluded because ²⁰Ne comes from the initial
composition and takes no part in the burning. The pool holds X = 1.00e-02:

| nuclide | X | share |
|---|---|---|
| ¹⁸Ne | 5.407e-03 | 56.3% |
| ¹⁵O | 2.629e-03 | 27.4% |
| ¹⁴O | 1.568e-03 | 16.3% |

¹⁸Ne has a 2.4 s mean lifetime, and its decay is exactly what drives the ¹⁵O
side flow above. Restricted to the four nuclei of the diagnostic ratio, the
split is 63% ¹⁵O / 37% ¹⁴O, with the nitrogen isotopes negligible.

## 5. Sensitivity to the expansion timescale

| $\tau$ (s) | $R^{\rm final}$ | $f_{\rm enh}$ | $t_{\rm fo}$ (s) |
|---|---|---|---|
| 0.05 | 0.584 | 159 | 0.052 |
| 0.10 | 0.609 | 166 | 0.112 |
| 0.20 | 0.460 | 125 | 0.226 |
| 0.50 | 0.330 | 90 | 0.454 |
| 1.00 | 0.315 | 86 | 0.352 |

The final ratio is not monotonic in $\tau$ (`fig17`). It peaks near
$\tau \approx 0.1$ s and falls away on both sides. Fast cooling freezes the
composition before much $^{15}$O can be made; slow cooling leaves enough time
at moderate temperature for $^{15}{\rm N}({\rm p},\alpha)$ to destroy some of
what was made. The largest enhancement comes from the intermediate case, and
the whole range is only a factor of two wide — much smaller than the factor of
seven between the exponential and trajectory models.

## 6. Matched controls: what actually causes the factor of four

The two reference histories differ in **five** ways at once — the temperature
maximum, the density at that maximum, presence of a heating phase, time spent
hot, and the temperature–density relation. The τ-series in section 5 varies only
the cooling rate inside the `T9_0 = 0.20` family, so it cannot say which of the
five matters. A second series was run for that: exponential models pinned to the
trajectory's conditions **at its temperature maximum** (`T9_0 = 0.418`,
`rho_0 = 4.00e3`), with the expansion timescale as the only free variable.
Note `4.00e3` is the density *when the temperature peaks*, not the trajectory's
highest density — that is `2.200e4`, at the start.

| τ (s) | t above T₉=0.2 (s) | R final | f_enh |
|---|---|---|---|
| 0.20 | 0.44 | 0.437 | 119 |
| 0.70 | 1.55 | 1.182 | 321 |
| 2.00 | 4.42 | 1.998 | 544 |
| 7.00 | 15.48 | 1.959 | 533 |
| **17.35** | **38.37** | **1.603** | **436** |
| 20.00 | 44.23 | 1.527 | 415 |
| **trajectory** | **38.37** | **1.958** | **533** |

Plus one further single-variable run, `exp_peakT_only`, which raises the peak
temperature to 0.418 while keeping the reference density and τ.

**The τ = 17.35 s row is constructed, not sampled — and it matters.** For the
exponential law the time above a threshold is `3 τ ln(T9_0 / T9_t)`, so the
timescale reproducing the trajectory's own exposure can be solved for exactly:
τ = 17.352 s gives 38.37 s against the trajectory's 38.37 s. Without it the
decomposition would have to interpolate across the τ = 7 → 20 s gap, a factor
of three in exposure.

It also defuses a trap. `exp_matched_tau7p0` returns **1.959** against the
trajectory's **1.958** — agreement to 0.04%. Read naively that says the
exponential family reproduces the trajectory exactly and the residual is 1.00.
It does not: that run reaches the same value at **less than half** the
trajectory's exposure. R is non-monotonic in exposure, peaking near τ = 2 s and
declining after, so the family crosses the trajectory's value twice for reasons
unconnected with matching. The comparison has to be made at equal exposure.

**A sequential factorisation** — each row changes exactly one thing, *in the
order listed*:

| Change | R before | R after | Factor |
|---|---|---|---|
| Peak temperature, 0.200 → 0.418 | 0.460 | 0.638 | ×1.39 |
| Density at T₉ max, 1.5e4 → 4.00e3 | 0.638 | 0.437 | ×0.69 |
| Exposure, 0.44 → 38.37 s | 0.437 | 1.603 | ×3.66 |
| Heating phase + T–ρ path | 1.603 | 1.958 | ×1.22 |
| **Combined** | **0.460** | **1.958** | **×4.26** |

**The individual factors are order-dependent.** The network is non-linear, so
the temperature and density steps interact. Running them the other way round
(`exp_rho_only`: density first, at the reference temperature):

| Order | temperature step | density step |
|---|---|---|
| T first, then ρ | ×1.39 | ×0.69 |
| ρ first, then T | ×0.73 | ×1.30 |

Both paths share their endpoints and their combined ×0.951 — they must. But the
individual factors **reverse direction**: lowering the density *raises* R at
T₉ = 0.20 and *lowers* it at T₉ = 0.418. These are not independent contributions
and cannot be recombined in another order. Exposure is the largest term under
both orderings tested — note it is the final step in both, so this is not a test
of every possible ordering.

Three further things follow.

**Peak temperature and density very nearly cancel.** Their combined factor is
×0.951 — *below* unity. Matching the trajectory's peak conditions alone leaves
the ratio slightly under the reference exponential model's. Essentially the
entire factor of 4.26 now rests on exposure (×3.66) and the residual (×1.22).

**Exposure alone does not explain it.** An exponential model matching the
trajectory in temperature maximum, density at that maximum *and* exposure time
reaches only 1.603 against the trajectory's 1.958. The residual ×1.22 belongs to
the two things an exponential prescription cannot reproduce: the heating phase,
during which the composition is already partly processed before maximum
temperature, and the trajectory's own T–ρ path, which is not the `T ∝ ρ^(1/3)`
of the analytic law. It is much smaller than the exposure term, but it has not
vanished.

**The exposure dependence is non-monotonic**, peaking near 4–5 s. Past that,
`15N(p,α)12C` keeps destroying A=15 material after production has stopped — the
same mechanism as in the cooler τ-series, at a different scale.

See `figures/fig18_matched_control.png`.

## 7. Network size

| Network | Range | Nuclides | Reactions | $R^{\rm final}$ | $\Delta R$ |
|---|---|---|---|---|---|
| Small | $Z\leq10$ | 68 | 535 | 1.9584 | +0.162% |
| Intermediate | $Z\leq20$ | 201 | 2119 | 1.9552 | −0.001% |
| Large | $Z\leq30$ | 370 | 4172 | 1.9552 | — |

**The diagnostic ratio is controlled almost entirely by local CNO cycling.**
Going from 68 nuclides to 370 changes it by 0.16%, so leakage into heavier
reaction cycles does not affect the CNO isotopic outcome at the level that
matters here, and the small network remains adequate for this diagnostic.

**The 0.16% is far above numerical noise.** Repeating `traj_ref` at `rtol`
1e-6 and 1e-7 gives 1.9584001529 and 1.9584000382 against the production
1.9584000261, so the solver reproduces $R^{\rm final}$ to 6.5e-8 relative. The
z10-to-z30 difference (1.6e-3) is 25000 times that, and even the z20-to-z30
difference (1.2e-5) is 185 times it, so the sign of the latter is resolved.

**The convergence is not uniform.** The two larger networks agree with each
other to 0.001%, so essentially the whole 0.16% sits between 68 and 201
nuclides.

**The cause is not what it looks like.** It is *not* the CNO-II/III path of
section 4a: every nuclide in that path has Z ≤ 10 and is present in all three
networks. Nor is it direct leakage — the largest reaction flow leaving the
`Z ≤ 10` network during the burning episode, integrated over t = 60–200 s, is
6 orders of magnitude below the cycle throughput (3.9e-6 of it).

It is an indirect effect through hydrogen. The `Z ≤ 10` network has no
proton-capture product for ²⁰Ne, so ²⁰Ne survives at X = 1.4e-3 instead of
being burnt to 6.3e-6 as it is in `Z ≤ 20`. The Ne–Na and Mg–Al processing that
follows in the larger networks consumes an extra 3.8e-3 in hydrogen mass
fraction — 0.53% of the total — building 5.85e-3 of material above neon, of
which the small network has exactly none. The truncated network therefore burns
at a slightly higher proton abundance, lifting every CNO abundance: A=15 by
0.42%, A=14 by 0.25%. The two groups do not respond equally, so their ratio
comes out 0.16% higher.

That does not mean the heavier nuclei do nothing — and on this profile they do
a great deal, because the material spends 38 s above $T_9 = 0.2$. Final-to-initial mass-fraction ratios from
the $Z\leq30$ run:

| element | initial | final | factor |
|---|---|---|---|
| Ne | 1.55e-3 | 6.29e-6 | ×0.004 |
| Na | 3.40e-5 | 1.59e-5 | ×0.47 |
| Mg | 6.50e-4 | 3.62e-6 | ×0.006 |
| Al | 5.60e-5 | 1.24e-6 | ×0.02 |
| Si | 6.70e-4 | 2.37e-5 | ×0.035 |
| S | 3.40e-4 | 1.16e-4 | ×0.34 |
| Ar | 8.70e-5 | 1.41e-4 | ×1.62 |
| Ca | 6.40e-5 | 4.71e-3 | **×73.6** |
| Fe | 1.30e-3 | 1.24e-3 | ×0.95 |

This is much stronger than the Ne-Na and Mg-Al cycling the proposal anticipates:
neon, magnesium, aluminium and silicon are almost entirely consumed and the
material piles up on calcium, which gains 4.6e-3 in mass fraction — comparable
to the whole CNO budget. The losses from Ne through S, roughly 3-4e-3, account
for most of that gain, so the flow is a one-way run up the chain rather than
cycling. Iron is essentially untouched (×0.95), as it must be at these
temperatures.

**Treat the calcium result with caution — but not for the reason of network
truncation.** Calcium sits at Z = 20, well inside the $Z\leq30$ network; the
abundance falls away above it (Sc 1.7e-4, Ti 5.6e-4, V 1.5e-5, Cr 8.7e-6) and
zinc at the network boundary is not populated at all, so the network edge is
demonstrably not binding. What does warrant caution is that this is a single
zone at fixed composition with no convective mixing, burnt for 38 s above
$T_9 = 0.2$ — an exposure that is itself an analytic assumption — and that
convergence against a network larger than $Z\leq30$ has not been tested. It is
reported because it is what the calculation gives, not as a prediction for nova
ejecta.

None of it feeds back on the $A=14$ and $A=15$ CNO isotopes at the level that
matters: extending the network changes $R_{15/14}$ by 0.16%, and almost all of
that comes from the CNO-II/III truncation at $Z\leq10$, not from the heavy-element
flow.

Practically, this means the small network is the one to use for the remaining
work: it gives the same CNO answer to 0.16% roughly six times faster (877 s against
5731 s for one trajectory run). The larger networks are still worth running
when the question is about Ne-Na or Mg-Al, where they are the whole point.

## 8. Numerical quality

| Check | Result |
|---|---|
| Baryon number, every reaction | exact, all three networks |
| Charge changes | only $\pm1$, i.e. beta decays with the leptons left out of the ReacLib record |
| Dead-end nuclides | none |
| Infinite or prompt rates in range | none |
| Beta-decay rates vs measured half-lives | within 0.2% for $^{13}$N, $^{14}$O, $^{15}$O, $^{17}$F, $^{18}$F |
| NucNetPy XML export vs the archive | identical counts, rates and masses |
| $\sum_i A_i Y_i$ at the end of a run | 1 to within 1.7e-14, every run |
| All seventeen runs | solver reported success and reached $t_{\rm end}$ |
| Network-size convergence | 0.16% between the smallest and largest network |
| Tolerance convergence | measured: $R^{\rm final}$ = 1.9584001529 / 1.9584000382 / 1.9584000261 at `rtol` 1e-6 / 1e-7 / 1e-8, i.e. reproducible to 6.5e-8 relative |

Four repairs were needed to the network before it behaved:

1. **Particle-unbound nuclides were dead ends.** ReacLib produces $^{5}$Li,
   $^{8}$Be and $^{9}$B but gives them no way to decay, so they slowly
   accumulated material that should have gone back to $\alpha$ particles and
   protons. Every reaction producing one is now written directly to its
   break-up products, the way $^{8}$B beta decay is normally written as
   $^{8}{\rm B} \rightarrow 2\alpha$.

2. **One prompt proton emitter was left in.** $^{18}{\rm Na} \rightarrow
   {\rm p} + ^{17}{\rm Ne}$ carries a rate of $5\times10^{20}$ s$^{-1}$ in
   ReacLib. It was eliminated the same way. This is found automatically now,
   by testing for one-body rates above $10^{12}$ s$^{-1}$.

3. **The ReacLib fits must not be extrapolated below $T_9 = 0.01$.** Below it
   some neutron-induced rates run away to $10^{24}$ s$^{-1}$ and beyond. Both
   thermodynamic histories are floored at $T_9 = 0.01$, which costs nothing:
   charged-particle reactions are already dead there, and the beta decays that
   still matter during the cooling tail do not depend on temperature.

4. **The neutron-richness cut was throwing away stable isotopes.** A flat
   $N - Z \leq 3$ band is fine for the CNO region but excludes $^{56}$Fe,
   $^{44}$Ca and $^{36}$S, because the valley of stability bends away from
   $N = Z$ as charge rises. Those isotopes are in the starting composition, so
   their mass was being silently redistributed over the survivors. The cut is
   now $N - Z \leq 3 + Z/4$, which retains every stable isotope of every
   element in the composition, and what a network cannot hold is reported
   explicitly in each run's `composition_truncation` record — 0.32% of the mass
   for $Z\leq10$, 0.13% for $Z\leq20$, nothing for $Z\leq30$.

None of the four changed the diagnostic ratio. $R^{\rm final}$ for the
exponential model was 0.46004 before any of them and 0.46004 after all of them.
That is reassuring rather than surprising: renormalising a composition scales
every abundance by the same factor, and unbound light nuclei and neutron-rich
iron do not talk to the CNO cycle. Each was still a defect that would have been
fair to ask about.

## 9. What is a modelling assumption, not a result

- **The trajectory is an input; the integration span is a choice.** The
  profile `data/trajectories/iliadis2002_S1_synthetic_benchmark.txt` is used as
  given, and it is a synthetic construction rather than hydrodynamic output —
  see the note at the top of this file. The runs stop at its last row,
  $t = 3000$ s, rather than the $3.15\times10^7$ s named in the proposal,
  because continuing would mean holding the final temperature and density
  fixed for four further decades of time. That is an extrapolation, not a
  result, and nothing turns on it: the ratio is within one per cent of its
  final value by 144.5 s.
- **The temperature is floored at $T_9 = 0.01$** because ReacLib's fits are not
  made for lower temperatures. For this profile the floor never binds: the
  history is generated with the same floor and approaches it asymptotically
  from above, so no tabulated value is altered. Only beta decays are still
  running by then, and they do not depend on temperature.
- **Deuterium is set to zero** in the initial composition. Material accreted
  from the companion has no primordial deuterium left; including the solar
  value puts a spurious spike on the energy generation in the first
  microsecond and changes nothing else.
- **No screening.** At $T_9 \leq 0.45$ and $\rho \leq 2.2\times10^4$ g
  cm$^{-3}$ the plasma is weakly coupled and screening changes the CNO rates by
  well under one per cent. NucNetPy's `SkyNetScreening` can be switched on if a
  quantitative bound is wanted.
- **Solar starting composition**, as the proposal specifies. Real nova envelopes
  are mixed with white-dwarf material, which would raise the CNO abundances by
  a large factor and change the energetics, though not obviously the ratio.
