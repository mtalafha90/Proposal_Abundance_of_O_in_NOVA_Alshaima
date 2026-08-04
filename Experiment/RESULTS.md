# Results

Every number below comes from the runs in `results/`, computed with NucNetPy.
Figures are in `figures/`; the LaTeX versions of the tables are in
`results/tables/`.

The trajectory model uses the measured nova profile
`data/trajectories/nova_profile_rescaled.txt`, which starts at
$T_9 = 0.09128$, $\rho = 2.211\times10^4$ g cm$^{-3}$, peaks at
$T_9 = 0.4481$ at $t = 103.89$ s, and is followed to the end of the file at
$t = 1.13\times10^5$ s.

## 1. The two reference calculations

| Model | Thermodynamic behaviour | $R^{\rm initial}$ | $R^{\rm final}$ | $f_{\rm enh}$ | $t_{\rm fo}$ | $T_{9,\max}$ |
|---|---|---|---|---|---|---|
| Exponential | immediate cooling | 3.68e-3 | **0.460** | 125 | 0.23 s | 0.200 |
| Trajectory | delayed temperature peak | 3.68e-3 | **3.230** | 879 | 155 s | 0.448 |

Both agree with the values quoted in the proposal, which gives
$R^{\rm final} \approx 0.4$–$0.5$ for the exponential model and "of order a
few" for the trajectory model. The shapes of the curves also match the
proposal's figures: an early bump, a sharp dip, a long climb, an overshoot, and
a flat freeze-out plateau (`figures/fig07_ratio_comparison.png`).

The ratio peaks at 3.97 during the burning phase and settles at 3.23, so the
trajectory model enhances $^{15}$N-rich material by a factor of 879 over the
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
| $^{14}$N(p,γ) — feeds A=15 | 125.0 s (T₉ 0.178) | 0.14 s |
| $^{12}$C(p,γ) — feeds the cycle | 128.6 s (T₉ 0.162) | 0.23 s |
| $^{15}$N(p,α) — drains A=15 | 160.0 s (T₉ 0.083) | 0.59 s |
| **numerical $t_{\rm fo}$** | **155.3 s** | **0.23 s** |

In both models $t_{\rm fo}$ falls inside the interval the crossings span, close
to the last one — the reaction that keeps changing the ratio longest. That is
the honest relationship: the crossing interprets $t_{\rm fo}$, it does not
define it.

**How the thermodynamic timescales are computed.** $T_9$ and $\rho$ come from
the same thermo callable the solver uses, sampled at the output times; for the
trajectory that is linear interpolation between table rows. No smoothing,
anywhere. Derivatives are `numpy.gradient` — second-order central differences on
the non-uniform grid, second-order one-sided at the ends. The table is finer
than the output grid through the burning episode (0.018 s vs 0.043 s), so
differentiating the interpolant is well-posed. $\tau_T$ diverges where
d$T_9$/d$t$ = 0 (at the peak and at minor stationary points); it is not
regularised, is drawn with a ceiling in the figures, and the nuclear-versus-
thermodynamic comparison is only used on the cooling branch where d$T_9$/d$t$ is
large — which is where freeze-out happens anyway.

**The trajectory produces more $^{15}$N because it burns hotter, for longer.**
The peak temperature is 0.448 against 0.200, and the measured profile holds the
material above $T_9 = 0.2$ for 17.7 seconds and above $T_9 = 0.1$ for 75 seconds
rather than cooling away from the peak immediately. Those are two of the five
things that differ between the histories; section 6 separates them with matched
controls and finds the exposure time to be the largest single contribution.

## 3. What happens to oxygen

The proposal is about oxygen in novae, so it is worth recording the elemental
outcome as well as the isotopic ratio. Total elemental mass fractions, start to
finish:

| Element | Initial | Exponential | Trajectory ($Z\leq10$) | Trajectory ($Z\leq30$) |
|---|---|---|---|---|
| C | 2.30e-3 | 1.11e-4 (×0.05) | 1.04e-3 (×0.45) | 1.04e-3 (×0.45) |
| N | 8.00e-4 | 2.27e-3 (×2.8) | 7.77e-3 (×9.7) | 7.74e-3 (×9.7) |
| O | 5.79e-3 | 6.78e-3 (×1.2) | 4.28e-6 (×0.0007) | 1.24e-6 (×0.0002) |

**The two thermodynamic histories do opposite things to oxygen.** The
exponential model, which never gets above $T_9 = 0.2$, leaves oxygen slightly
enhanced: it burns carbon into nitrogen and a little oxygen, and then freezes.
The trajectory model, which reaches $T_9 = 0.448$, destroys oxygen by a factor
of 1350 and converts it into nitrogen, which ends up ten times its solar value.
$^{16}$O falls from a mass fraction of 5.77e-3 to 2.96e-6.

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

Five of the six carry essentially no traffic outside the cycle. **¹⁵O does.**
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

**The limiting step changes with temperature**, so the classical CNO
bottleneck does not hold throughout. The proton-capture timescale of
$^{14}$N is 244 s at the starting conditions, longer than the beta lifetimes
of $^{14}$O (101.9 s) and $^{15}$O (176.0 s), so $^{14}$N(p,$\gamma$) limits the
flow during the quiescent phase. It drops below the $^{15}$O lifetime at
$T_9 \simeq 0.10$ (t $\simeq$ 74 s) and reaches $3.9\times10^{-5}$ s at the
temperature maximum — six orders of magnitude faster. Through the whole burning
episode the cycle is therefore **beta-limited**.

**There are three waiting points, and ¹⁵O is not the largest.** At the
temperature maximum (t = 103.907 s, T₉ = 0.4467), 96% of the circulating
material — C, N, O and F isotopes plus ¹⁸Ne and ¹⁹Ne — sits on three
beta-unstable nuclei:

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

## 6. Matched controls: what actually causes the factor of seven

The two reference histories differ in **five** ways at once — the temperature
maximum, the density at that maximum, presence of a heating phase, time spent
hot, and the temperature–density relation. The τ-series in section 5 varies only the cooling
rate inside the `T9_0 = 0.20` family, so it cannot say which of the five
matters. A second series was run for that: exponential models pinned to the
trajectory's conditions **at its temperature maximum** (`T9_0 = 0.4481`,
`rho_0 = 4.07e3`), with the expansion timescale as the only free variable.
Note `4.07e3` is the density *when the temperature peaks*, not the trajectory's
highest density — that is `2.211e4`, near the start.

| τ (s) | t above T₉=0.2 (s) | R final | f_enh |
|---|---|---|---|
| 0.2 | 0.48 | 0.537 | 146 |
| 0.7 | 1.69 | 1.343 | 365 |
| 2.0 | 4.84 | 2.009 | 547 |
| 7.0 | 16.94 | 1.951 | 531 |
| 20.0 | 48.40 | 1.524 | 415 |
| **trajectory** | **17.68** | **3.230** | **879** |

Plus one further single-variable run, `exp_peakT_only`, which raises the peak
temperature to 0.4481 while keeping the reference density and τ.

**A sequential factorisation** — each row changes exactly one thing, *in the
order listed*:

| Change | R before | R after | Factor |
|---|---|---|---|
| Peak temperature, 0.200 → 0.448 | 0.460 | 0.697 | ×1.52 |
| Density at T₉ max, 1.5e4 → 4.07e3 | 0.697 | 0.537 | ×0.77 |
| Exposure, 0.48 → 16.9 s | 0.537 | 1.951 | ×3.64 |
| Heating phase + T–ρ path | 1.951 | 3.230 | ×1.66 |
| **Combined** | **0.460** | **3.230** | **×7.02** |

**The individual factors are order-dependent, and severely so.** The network is
non-linear, so the temperature and density steps interact. Running them the
other way round (`exp_rho_only`: density first, at the reference temperature):

| Order | temperature step | density step |
|---|---|---|
| T first, then ρ | ×1.52 | ×0.77 |
| ρ first, then T | ×0.89 | ×1.31 |

Both paths share their endpoints and their combined ×1.166 — they must. But the
individual factors don't merely differ in size, they **reverse direction**:
lowering the density *raises* R at T₉ = 0.20 and *lowers* it at T₉ = 0.448.
So these are not independent contributions and cannot be recombined in another
order. The exposure term is the largest under both orderings tested — note it
is the final step in both, so this is not a test of every possible ordering —
which is why the headline conclusion survives.

Three further things follow.

**Exposure time is the largest single contribution** (×3.64), ahead of peak
temperature (×1.52). The paper's thesis survives — but it needed this to be
shown rather than asserted.

**Exposure alone does not explain it.** An exponential model matching the
trajectory in temperature maximum, density at that maximum *and* exposure time reaches only
1.951 against the trajectory's 3.230. The residual ×1.66 belongs to the two
things an exponential prescription cannot reproduce: the heating phase, during
which the composition is already partly processed before maximum temperature,
and the trajectory's own T–ρ path, which is not the `T ∝ ρ^(1/3)` of the
analytic law.

**The exposure dependence is non-monotonic here too**, peaking near 5 s. Past
that, `15N(p,α)12C` keeps destroying A=15 material after production has
stopped — the same mechanism as in the cooler τ-series, at a different scale.

See `figures/fig18_matched_control.png`.

## 7. Network size

| Network | Range | Nuclides | Reactions | $R^{\rm final}$ | $\Delta R$ |
|---|---|---|---|---|---|
| Small | $Z\leq10$ | 68 | 535 | 3.2299 | +0.016% |
| Intermediate | $Z\leq20$ | 201 | 2119 | 3.2297 | +0.008% |
| Large | $Z\leq30$ | 370 | 4172 | 3.2294 | — |

**The diagnostic ratio is controlled entirely by local CNO cycling.** Going
from 68 nuclides to 370 changes it by less than two parts in ten thousand.
This settles the question the proposal poses: leakage into heavier reaction
cycles does not affect the CNO isotopic outcome, and the small network is
enough for this diagnostic.

That does not mean the heavier nuclei do nothing. The larger networks show
real processing above neon during the burning phase — neon down to a tenth of
its starting value, sodium up fourteenfold, sulphur up threefold, argon up
sevenfold — which is the Ne-Na and Mg-Al cycling the proposal anticipates. Iron
is untouched, as it must be at these temperatures. None of it feeds back on the
$A=14$ and $A=15$ CNO isotopes: the flow runs one way, out of the CNO region and
upwards, and it is far too small to affect the CNO abundances it leaves behind.

Practically, this means the small network is the one to use for the remaining
work: it gives the same CNO answer roughly six times faster (1280 s against
7545 s for one trajectory run). The larger networks are still worth running
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
| $\sum_i A_i Y_i$ at the end of a run | 1 to within 2e-14, every run |
| All nine runs | solver reported success and reached $t_{\rm end}$ |
| Network-size convergence | 0.016% between the smallest and largest network |
| Tolerance convergence | $R^{\rm final}$ stable to six figures between `rtol` 1e-6 and 1e-8 |

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

- **The trajectory is measured data; the integration span is a choice.** The
  profile `data/trajectories/nova_profile_rescaled.txt` is used as given. The
  runs stop at its last row, $t = 1.13\times10^5$ s, rather than the
  $3.15\times10^7$ s named in the proposal, because continuing would mean
  holding the final temperature and density fixed for two further decades of
  time. That is an extrapolation, not a result, and nothing turns on it: the
  ratio is within one per cent of its final value by 155 s and the last
  relevant decay, $^{18}$F, is complete by $10^5$ s.
- **The temperature is floored at $T_9 = 0.01$**, which the profile reaches at
  $t = 344$ s, because ReacLib's fits are not made for lower temperatures. Only
  beta decays are still running by then, and they do not depend on temperature.
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
