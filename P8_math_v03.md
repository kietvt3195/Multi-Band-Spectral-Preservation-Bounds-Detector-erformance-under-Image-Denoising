# P8 Math v0.3 — Corrected Theoretical Core (Theorem 1, sound proof)

**Status**: 16/05/2026 | Supersedes `P8_math_v02.md` (15/05/2026)
**Author**: Võ Thành Kiệt — VSB-TUO | **Target venue**: IEEE TIP
**Scope of this revision**: theory only. Empirical section is a stub pending
Notebook 1 v4 + Notebook 3 outputs (see §8). Numbers from v02 that came from
the stale 5-image run are NOT carried over.

---

## §0. WHY v0.3 EXISTS — what was wrong in v0.2

A line-by-line audit of `bfpi_lib_v02.py`, `P8_math_v02.md`, the paper
outline, and `paper8/sections/{03_theory,99_appendix}.tex` found **three
distinct defects** in the v0.2 theoretical core. v0.3 fixes all three.

### D1 — Two incompatible definitions of the theorem quantity

Two "camps" defined the symmetric preservation differently:

- **Code camp** (`bfpi_lib_v02.py`, `P8_math_v02.md` Def 1′):
  $\rho_k = 1 - \min(1,\,|E_k(\hat y)-E_k(y)|/E_k(y))$
- **Paper camp** (outline §3.2, `03_theory.tex` `def:rho_sim`, appendix):
  $\rho_k = \min(\rho^{ratio}_k,\,1/\rho^{ratio}_k)$

These agree on the over-smoothing branch ($\rho^{ratio}\le 1$) but **diverge on
the artifact branch** ($\rho^{ratio}>1$): e.g. at ratio $=2$ the code gives
$0$, the paper gives $0.5$. The artifact branch is exactly where DnCNN lives
under domain shift, so this is not a corner case.

### D2 — The per-band inequality is FALSE (the central error)

Both camps' $\rho_k$ are functions of band **energies** only — i.e. of the
**magnitude spectrum** $|\hat x|^2$. They discard Fourier phase. But the proof
chain runs through $\|y-\hat y\|_2^2$, which via Parseval equals
$\tfrac1{HW}\|\hat Y-\hat Y'\|_2^2$ — a quantity that depends on the **full
complex spectrum, phase included**.

A phase-blind quantity cannot upper-bound a phase-sensitive one.
**Counterexample**: let $\hat y$ have the same band energies as $y$ in every
band but opposite phase. Then $\rho^{ratio}_k=1$, so $\rho_k=1$ and the v0.2
right-hand side $E_k(y)(1-\rho_k)^2 = 0$, while the true per-band distance
$\|\hat Y_k-\hat Y'_k\|^2$ is large. The v0.2 `Lemma perband` is therefore not
"uncertain" — it is **not true** without an unstated phase-alignment
assumption. (Even granting that assumption, its final algebraic step
$|1-\sqrt{\rho^{ratio}}|\le 1-\rho_k$ still fails for $\rho^{ratio}>1$, since
it reduces to $s+1/s\le 2$ with $s=\sqrt{\rho^{ratio}}>1$.)

Consequence: the v0.2 "0/25 violations" result is a **post-hoc calibration of
$C$**, not a verification of the proof. Fitting $C$ to the data absorbs all
slack, including the dropped phase term.

### D3 — The K-band partition does not cover the frequency domain

`get_band_boundaries` uses `linspace(0, π, K+1)`, so the bands tile only the
disk $R\le\pi$. The 2-D DFT grid is the square $[-\pi,\pi]^2$, on which the
radial magnitude $R=\sqrt{u^2+v^2}$ reaches $\sqrt2\,\pi\approx 4.44$ at the
corners. Frequencies with $R\in(\pi,\sqrt2\,\pi]$ fall in **no band**, so the
band-decomposition step $\|\hat Y-\hat Y'\|^2=\sum_k\|\cdot\|_{B_k}^2$ is not
an equality — it omits the corner energy. Fix: the last band is **unbounded**.

---

## §1. CORRECTED DEFINITIONS

Notation: $y$ = clean image, $\hat y$ = denoised image, both in
$\mathbb R^{H\times W}$ on the luminance (Y) channel; $\hat Y=\mathrm{DFT}(y)$,
$\hat Y'=\mathrm{DFT}(\hat y)$. Restriction to band $B_k$ is written
$\hat Y_k:=\hat Y|_{B_k}$.

### Definition 1 — K-band radial partition (last band unbounded)

For spatial frequency $(u,v)$, let $R(u,v)=\sqrt{u^2+v^2}$. The **K-band
partition** is
$$
B_k := \{(u,v): R\in[(k-1)\pi/K,\ k\pi/K)\},\quad k=1,\dots,K-1,
$$
$$
\boxed{\,B_K := \{(u,v): R\ge (K-1)\pi/K\}\,}
$$
i.e. the **last band is unbounded** and absorbs the corner frequencies. With
this, $\{B_k\}_{k=1}^K$ is a genuine partition of the entire DFT grid. We fix
$K=4$, giving bands $[0,\tfrac\pi4),[\tfrac\pi4,\tfrac\pi2),[\tfrac\pi2,\tfrac{3\pi}4),[\tfrac{3\pi}4,\infty)$.

> **Change vs v0.2**: band 4 was $[\tfrac{3\pi}4,\pi)$; it is now
> $[\tfrac{3\pi}4,\infty)$. Required for Lemma B (§3) to be an equality.

### Definition 2 — Band energy

$$E_k(x):=\sum_{(u,v)\in B_k}|\hat x(u,v)|^2=\|\hat X_k\|_2^2,\qquad
  E(x)=\sum_{k=1}^K E_k(x).$$

### Definition 3 — Band-wise spectral preservation $\rho_k$ (phase-aware) ★

This is the quantity that enters Theorem 1. Define the **normalized band
distortion**
$$
\delta_k(y,\hat y):=\frac{\|\hat Y_k-\hat Y'_k\|_2}{\|\hat Y_k\|_2}
  =\sqrt{\frac{\sum_{(u,v)\in B_k}|\hat Y(u,v)-\hat Y'(u,v)|^2}{E_k(y)}}\ \ge 0,
$$
and the **band-wise spectral preservation**
$$
\boxed{\ \rho_k(y,\hat y) := 1-\delta_k(y,\hat y)\ \le 1.\ }
$$

Properties:
- $\rho_k=1 \iff \hat Y_k=\hat Y_k$ exactly — perfect preservation of **both
  magnitude and phase** in band $k$.
- $\rho_k$ is **not** clamped: $\rho_k<0$ is possible and *meaningful* — it is
  the strong-distortion regime where the denoised band differs from the clean
  band by more than the clean band's own magnitude (heavy artifact injection
  or phase inversion).
- For plotting/description one may report $\rho_k^{+}:=\max(0,\rho_k)$, but the
  theorem and its proof use the **unclamped** $\rho_k$ (equivalently $\delta_k$).
- Defined when $E_k(y)>0$. For natural images this holds for all $k$; the
  degenerate case $E_k(y)=0$ is handled by taking that band's contribution as
  the raw distance $\|\hat Y_k-\hat Y'_k\|^2$ (see Lemma C remark).

> **This replaces** v0.2 Definition 1′ and `03_theory.tex` `def:rho_sim`. The
> old name "$\rho^{sim}$" / `\rhoSim` macro should be retired (see §7).

### Definition 4 — Energy ratio $\rho^{ratio}_k$ (diagnostic only)

$$\rho^{ratio}_k(y,\hat y):=\frac{E_k(\hat y)}{E_k(y)}\in[0,\infty).$$

Unchanged from `bfpi_lib_v02.py` `mode='ratio'`. **Not** used in Theorem 1.
Used only for two-mode failure classification (§6).

### Definition 5 — Band phase coherence $\gamma_k$, and the v0.2 reconciliation

$$\gamma_k(y,\hat y):=\frac{\mathrm{Re}\,\langle\hat Y'_k,\hat Y_k\rangle}
  {\|\hat Y'_k\|_2\,\|\hat Y_k\|_2}\in[-1,1].$$

Then, expanding the squared norm and writing $r=\rho^{ratio}_k$:
$$
\boxed{\ \delta_k^2 = 1 + r - 2\gamma_k\sqrt r
       = \underbrace{(1-\sqrt r)^2}_{\text{magnitude}}
       + \underbrace{2\sqrt r\,(1-\gamma_k)}_{\text{phase}}.\ }
$$

This identity makes the v0.2 error explicit: **v0.2's energy-only $\rho$
corresponds to assuming $\gamma_k\equiv 1$ (perfect phase preservation)** and
dropping the phase term $2\sqrt r(1-\gamma_k)$. v0.3's $\rho_k$ keeps the full
expression. $\gamma_k$ need not be computed for Theorem 1 (it is contained in
$\delta_k$); it is introduced only to (a) explain the correction and (b) enable
an optional phase-aware refinement of the two-mode diagnostic (§6, future work).

### Definition 6 — Detector Spectral Sensitivity Profile (DSSP) $\alpha_k$

For detector $D_\theta$ and band index $k$, the **band-$k$-masked** image is
$y^{(-k)}:=\mathcal F^{-1}(\hat Y\cdot\mathbb 1_{B_k^{\,c}})$. The unnormalized
sensitivity over a probe set $\{y_i\}_{i=1}^N$ is
$\tilde\alpha_k=\tfrac1N\sum_i[F_1(D_\theta(y_i))-F_1(D_\theta(y_i^{(-k)}))]$,
normalized to the simplex $\alpha_k=\tilde\alpha_k/\sum_j\tilde\alpha_j$,
$\sum_k\alpha_k=1$.

Unchanged from v0.2 Definition 2 / Notebook 2. **The masking in the $\alpha$
probe must use the same Definition-1 band edges** (including the unbounded
band $K$) as the $\rho_k$ computation — otherwise $\alpha_k$ and $\rho_k$
refer to different bands. (Notebook 2 currently uses `linspace(0,π,5)`; it
needs the same Definition-1 fix as Notebook 1 — see §8.)

### Definition 7 — Active denoising regime (empirical scoping)

A denoiser is **active** at noise level $\sigma$ if its mean band distortion
exceeds a threshold, $\overline{\sum_k\alpha_k\delta_k^2}>\tau_{\text{active}}$.
Near-identity transforms (e.g. a weak fixed Gaussian filter) fall below
$\tau_{\text{active}}$ and produce $\Delta F_1\approx0$ by training-augmentation
tolerance, independent of the bound. $\tau_{\text{active}}$ is to be calibrated
on the v4 data; the specific value in v0.2 is not carried over.

---

## §2. THEOREM 1 (corrected) — Spectral Preservation Bound

**Assumptions (stated explicitly, to be disclosed in the paper).**

- **(A1) Lipschitz detection score.** The map $y\mapsto F_1(D_\theta(y))$ is
  Lipschitz with constant $L$ on the relevant image domain:
  $|F_1(D_\theta(y))-F_1(D_\theta(\hat y))|\le L\,\|y-\hat y\|_2$.
  This is a *modeling assumption*: $F_1$ is not globally Lipschitz (it is
  built from thresholded, matching-based box counts). We assume locally
  Lipschitz behaviour on the noise/denoising perturbation manifold and
  **estimate $L$ empirically** (it is folded into $C$ below).
- **(A2) Single global constant.** $C$ as derived below is image-dependent;
  the paper reports one fitted $C$ across images, which is an approximation
  (Remark 1).

**Theorem 1.** Under (A1), for any clean image $y$ and denoised image $\hat y$,
$$
\boxed{\ \bigl|F_1(D_\theta(y))-F_1(D_\theta(\hat y))\bigr|^2
   \ \le\ C\sum_{k=1}^{K}\alpha_k(D_\theta)\,\bigl(1-\rho_k(y,\hat y)\bigr)^2\ }
$$
with $\rho_k$ the phase-aware preservation of Definition 3, $\alpha_k$ the
DSSP of Definition 6, and $C=\dfrac{L^2}{HW}\max_k\dfrac{E_k(y)}{\alpha_k}\ge0$.

> The boxed equation is **textually identical** to v0.2 / `03_theory.tex`. Only
> the *definition* of $\rho_k$ changed (Definition 3). This keeps the paper's
> Corollary 1, $\alpha_k$ machinery, and $K=4$ discussion intact.

---

## §3. PROOF OF THEOREM 1

**Lemma A (Parseval–Plancherel, 2-D DFT).** For $y,\hat y\in\mathbb R^{H\times W}$,
$\|y-\hat y\|_2^2=\tfrac1{HW}\|\hat Y-\hat Y'\|_2^2.$ *(Standard.)*

**Lemma B (Exact band decomposition).** If $\{B_k\}_{k=1}^K$ partitions the DFT
grid (Definition 1, last band unbounded), then
$$\|\hat Y-\hat Y'\|_2^2=\sum_{k=1}^K\|\hat Y_k-\hat Y'_k\|_2^2.$$
*Proof.* Each frequency lies in exactly one $B_k$, so the restrictions
$(\hat Y-\hat Y')|_{B_k}$ have pairwise disjoint support; the squared norm is
additive over disjoint supports. The **unbounded** $B_K$ is what makes the
$B_k$ cover every grid frequency, including the corners $R\in(\pi,\sqrt2\,\pi]$;
with the v0.2 bounded band 4 this step omitted the corner energy. $\square$

**Lemma C (Per-band identity).** With $\delta_k,\rho_k$ as in Definition 3 and
$E_k(y)>0$,
$$\|\hat Y_k-\hat Y'_k\|_2^2 = E_k(y)\,\delta_k^2 = E_k(y)\,(1-\rho_k)^2.$$
*Proof.* Immediate from $\delta_k=\|\hat Y_k-\hat Y'_k\|_2/\|\hat Y_k\|_2$ and
$\|\hat Y_k\|_2^2=E_k(y)$. This is an **equality by construction** — it
replaces the v0.2 `Lemma perband`, which was a false inequality. *Remark:* if
$E_k(y)=0$, take the band's contribution to be the raw
$\|\hat Y_k-\hat Y'_k\|_2^2$ and leave $\rho_k$ undefined for that band; the
main proof below still goes through term-by-term. $\square$

**Lemma D (Lipschitz of the F1 score).** This is Assumption (A1) restated;
it is *assumed*, not derived (see §2 and Remark 2).

**Main proof.**
$$
\begin{aligned}
|\Delta F_1|^2
&\le L^2\,\|y-\hat y\|_2^2 && \text{(Lemma D / A1)}\\
&= \frac{L^2}{HW}\,\|\hat Y-\hat Y'\|_2^2 && \text{(Lemma A)}\\
&= \frac{L^2}{HW}\sum_{k=1}^K\|\hat Y_k-\hat Y'_k\|_2^2 && \text{(Lemma B)}\\
&= \frac{L^2}{HW}\sum_{k=1}^K E_k(y)\,(1-\rho_k)^2 && \text{(Lemma C)}\\
&= \frac{L^2}{HW}\sum_{k=1}^K \alpha_k\cdot\frac{E_k(y)}{\alpha_k}\cdot(1-\rho_k)^2\\
&\le \frac{L^2}{HW}\,\Bigl(\max_k\frac{E_k(y)}{\alpha_k}\Bigr)
      \sum_{k=1}^K \alpha_k\,(1-\rho_k)^2 && (\ast)\\
&= C\sum_{k=1}^K \alpha_k\,(1-\rho_k)^2,\qquad
   C:=\frac{L^2}{HW}\max_k\frac{E_k(y)}{\alpha_k}.
\end{aligned}
$$
Step $(\ast)$ uses $\alpha_k>0$, $\sum_k\alpha_k=1$: a convex combination is
$\le$ its maximum term. $\square$

**Corollary 1 (Worst-case F1 loss).**
$$|\Delta F_1|\ \le\ \sqrt{\,C\sum_{k=1}^K\alpha_k(1-\rho_k)^2\,}.$$
A closed-form, detector-free safety guarantee: computable from a single image
pair $(y,\hat y)$ and a once-estimated $\alpha$ profile.

**Remark 1 (image-dependence of $C$).** $C$ depends on $y$ through
$\max_k E_k(y)/\alpha_k$ — band 1 (which contains the DC term) typically
dominates. Reporting a single fitted $C$ across a test set is an approximation;
honest options are (a) report $C$ as a high quantile / supremum over the
sample, or (b) make $C$ image-conditional. **Disclose this** (it is the
honest version of CLAUDE.md lock J3 "$C$ fitted post-hoc").

**Remark 2 (looseness of $(\ast)$).** The $\max_k$ step is loose. A Jensen-type
or per-band constant would tighten it. Keep this as stated future work
(matches the existing appendix TODO).

**Remark 3 (per-image vs aggregate).** Theorem 1 is a **per-image** statement.
The empirical pipeline averages $\rho_k$ over ~50 images per $(\text{method},\sigma)$
and pairs it with a dataset-level $\Delta F_1$. The validation is therefore an
aggregate-level consistency check, not a per-image test. **Disclose** this as a
limitation (it belongs in §5.7).

---

## §4. WHAT IS SOUND vs ASSUMED vs PENDING

| Element | Status |
|---|---|
| Definitions 1–6 | Sound, self-consistent |
| Lemma A (Parseval) | Sound (standard) |
| Lemma B (band decomposition) | Sound **iff** band $K$ unbounded (Def 1) |
| Lemma C (per-band identity) | Sound — equality by construction |
| Main proof chain | Sound given (A1) |
| (A1) $F_1\!\circ\!D$ Lipschitz | **Assumption** — disclose; $L$ folded into $C$ |
| Single global $C$ | **Approximation** — disclose (Remark 1) |
| Per-image vs aggregate validation | **Limitation** — disclose (Remark 3) |
| Empirical $C$, correlations, $\tau_{\text{active}}$ | **Pending** NB1 v4 + NB3 (§8) |

The phase defect (D2) and the partition defect (D3) are **fixed**, not merely
disclosed. The three remaining items are legitimate, disclosed modeling
choices — appropriate for a TIP paper, not errors.

---

## §5. THEORETICAL CONTRIBUTIONS (updated list)

1. **Definition 3** — phase-aware band-wise spectral preservation $\rho_k$
   (replaces the energy-only v0.2 quantity).
2. **Definition 5** — exact magnitude/phase decomposition
   $\delta_k^2=(1-\sqrt{\rho^{ratio}_k})^2+2\sqrt{\rho^{ratio}_k}(1-\gamma_k)$,
   which shows the v0.2 quantity is the $\gamma_k\equiv1$ special case.
3. **Definition 6** — DSSP $\alpha_k$ (unchanged).
4. **Theorem 1** — closed-form upper bound on $|\Delta F_1|^2$, now with a
   sound proof under explicit assumption (A1).
5. **Corollary 1** — detector-free worst-case F1-loss guarantee.
6. **Two-mode framework** — over-smoothing vs artifact injection via
   $\rho^{ratio}_k$ (§6).

---

## §6. TWO-MODE FAILURE DIAGNOSTIC (unchanged role)

Theorem 1's $\rho_k$ does not reveal failure *direction*. The diagnostic uses
$\rho^{ratio}_k$ (Definition 4): **faithful** $\rho^{ratio}_k\in[0.8,1.2]\ \forall k$;
**over-smoothing** $\exists k:\rho^{ratio}_k<0.5$; **artifact injection**
$\exists k:\rho^{ratio}_k>1.5$; **mixed** otherwise. This is purely
descriptive and is unaffected by the Theorem 1 correction.
*Optional future refinement:* a low $\gamma_k$ (Definition 5) flags phase
scrambling — a third failure axis the energy ratio cannot see.

---

## §7. CHANGES REQUIRED IN THE PAPER (`paper8/`)

- **`03_theory.tex`** — replace `def:rho_sim` with Definition 3 (phase-aware);
  add Definition 5 (magnitude/phase split) and the unbounded-band-$K$ note in
  `def:partition`; update `prop:rho_sim_properties` ($\rho_k\le1$, may be $<0$;
  $\rho_k=1\iff$ exact band match).
- **`99_appendix.tex`** — replace the broken `Lemma perband` with Lemma C
  (one-line identity); make `Lemma parseval` carry the unbounded partition;
  remove the `\TODO{Verify bất đẳng thức cuối}` (resolved). State (A1)
  explicitly as a hypothesis of the theorem.
- **`macros.tex`** — retire `\rhoSim`; the theorem quantity is now `\rhoB{k}`.
  Optionally add `\deltaB{k}`, `\gammaB{k}`.
- **§5.7 / §6** — add the disclosure paragraphs for Remarks 1 and 3.
- **Abstract / §1** — wording is safe: the boxed bound is unchanged, so the
  contributions list and abstract claims still hold.

All edits are surgical `str_replace` per Standards A1 — no full-file rewrite.

---

## §8. CHANGES REQUIRED IN THE NOTEBOOKS

The phase-aware $\rho_k$ **cannot be recovered from the v4 CSVs**: the v4
schema stores band *energies* and the energy ratio, but $\rho_k$ needs
$\|\hat Y_k-\hat Y'_k\|^2$ (a complex-spectrum quantity). One more notebook
revision is needed.

**Notebook 1 → v5** (small, surgical):
1. In the BFPI library, add a function computing the **complex per-band
   distance** $\|\hat Y_k-\hat Y'_k\|^2$ directly from the two DFTs (this is
   *simpler* than the existing `band_energy` loop), and derive $\rho_k$,
   $\delta_k$ per Definition 3.
2. In `get_band_boundaries`, make band $K$ **unbounded** (Definition 1): the
   last band's mask becomes `radial_map >= boundaries[K-1]` with no upper bound.
3. Add CSV columns `rho_phase_k` (the Definition-3 $\rho_k$). Keep the existing
   `rho_ratio_k` (diagnostic) and `E_band{k}_clean`. The old energy-based
   `rho_k` columns may stay, relabelled as descriptive only.

**Notebook 2** — apply the same unbounded-band-$K$ fix to
`mask_band_in_luminance` so $\alpha_k$ and $\rho_k$ share the band definition.
Re-running the $\alpha$ probe is cheap and deterministic.

**Notebook 3** — `verify_theorem1_for_detector` should read `rho_phase_k`
instead of the energy `rho_k`. (Recall the merge-key check passed: the
288-run CSV uses `denoise_method ∈ {autoencoder,bm3d,cae_pso,dncnn,
gaussian_filter,noisy}`, matching NB1.)

**Good news for the run in progress.** The current NB1 v4 sweep is **not
wasted**: `rho_ratio_k` (the two-mode diagnostic) and `E_band{k}_clean` remain
valid and are exactly what §6 needs. Only the *theorem* quantity changes, and
that is a one-function addition for the v5 run.

---

## §9. EMPIRICAL VALIDATION — STUB (pending data)

To be filled once NB1 v4 (running 16-May) and NB1 v5 + NB3 complete:

- Theorem 1 verification with the **phase-aware** $\rho_k$: violation count,
  fitted $C$ (with the BM3D $\sigma{=}1$ anomaly disclosure), bound tightness.
- Expectation: phase-aware $\delta_k$ is generally $\ge$ the magnitude-only
  distortion, so $(1-\rho_k)^2$ is larger and the fitted $C$ should be
  **smaller and more stable** than v0.2's $C=10.71$ — a sounder *and* tighter
  result. To be confirmed, not assumed.
- Cross-dataset BM3D stability, DSSP profile, two-mode classification:
  re-confirm on v4/v5 data. The v0.2 numbers (5-image run) are not carried
  over.

---

## §10. DELIVERABLES STATUS

| Item | Status |
|---|---|
| Corrected definitions (1–7) | ✅ this document |
| Theorem 1 + sound proof | ✅ this document |
| Reconciliation of the two $\rho$ definitions (D1) | ✅ Definition 3 supersedes both |
| Phase defect (D2) | ✅ fixed via phase-aware $\rho_k$ |
| Partition defect (D3) | ✅ fixed via unbounded band $K$ |
| Paper edit list | ✅ §7 |
| Notebook edit list (v5) | ✅ §8 |
| Empirical re-validation | ⏳ pending NB1 v4 → v5 → NB3 |
