"""
Within-family polygenic-score design — separating DIRECT genetic effects from
GENETIC NURTURE (indirect parental effects), on simulated trios.

Application companion for the PhD position
"Social Science Genetics" (VU Amsterdam & Amsterdam UMC, MSCA), whose central
question is: how much of the polygenic-score association with life outcomes is a
direct effect of the child's own genes, versus an indirect effect of the
parents' genes acting through the rearing environment ("nature of nurture")?

Method (Kong et al., Science 2018): for each child, split the parental alleles
into TRANSMITTED (the child's own genotype) and NON-TRANSMITTED. Non-transmitted
alleles are not in the child, so they cannot act directly — yet they still
predict the child's outcome through the environment parents create. Regressing
the outcome on transmitted (T) and non-transmitted (NT) polygenic scores:
    coef(NT) = genetic nurture (indirect)
    coef(T)  = direct + nurture
    direct   = coef(T) - coef(NT)
A naive population PGS regression conflates the two and overstates "nature".

This is exactly an instrumental / natural-experiment identification — Mendelian
segregation randomises which parental allele is transmitted — which is why an
economics/econometrics background transfers directly.

Synthetic data only. Deps: numpy, scipy, matplotlib.
Run:  python3 within_family_demo.py
"""
from __future__ import annotations
import json, warnings
import numpy as np
warnings.filterwarnings("ignore", category=RuntimeWarning)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

RNG = np.random.default_rng(3)
N = 4000        # trios
M = 200         # causal SNPs
BETA_DIRECT = 0.30
BETA_NURTURE = 0.15


def simulate_trios():
    """Return transmitted PGS (=child), non-transmitted PGS, and child outcome."""
    maf = RNG.uniform(0.1, 0.5, M)
    w = RNG.normal(0, 1, M)                       # per-SNP effect weights

    # Parental haplotypes (one allele each), per SNP, for father and mother.
    f1 = RNG.binomial(1, maf, (N, M)); f2 = RNG.binomial(1, maf, (N, M))
    m1 = RNG.binomial(1, maf, (N, M)); m2 = RNG.binomial(1, maf, (N, M))

    # Mendelian transmission: each parent passes one haplotype at random.
    pick_f = RNG.random((N, M)) < 0.5
    pick_m = RNG.random((N, M)) < 0.5
    t_f = np.where(pick_f, f1, f2); nt_f = np.where(pick_f, f2, f1)
    t_m = np.where(pick_m, m1, m2); nt_m = np.where(pick_m, m2, m1)

    child = t_f + t_m                              # transmitted dosage (child genotype)
    nontrans = nt_f + nt_m                         # non-transmitted dosage

    def pgs(dosage):
        s = dosage @ w
        return (s - s.mean()) / s.std()

    T = pgs(child)          # transmitted PGS (standardised), = child's own PGS
    NT = pgs(nontrans)      # non-transmitted PGS (standardised)

    # Outcome in raw effect units so the regression recovers the true parameters:
    # direct effect acts through the child's genotype (T); genetic nurture acts
    # through the parents' genotype, i.e. through BOTH transmitted and non-
    # transmitted alleles (T + NT). T and NT are independent by Mendelian
    # segregation, so:  y = (direct + nurture)*T + nurture*NT + noise.
    # We deliberately do NOT re-standardise y, so coef(T)=direct+nurture and
    # coef(NT)=nurture come out in the stated units.
    y = (BETA_DIRECT + BETA_NURTURE) * T + BETA_NURTURE * NT + RNG.normal(0, 1, N)
    return T, NT, y


def ols(X, y):
    """OLS with HC0-ish SE; returns betas, ses, p-values. X includes intercept col."""
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    sigma2 = (resid @ resid) / dof
    cov = sigma2 * XtX_inv
    se = np.sqrt(np.diag(cov))
    p = 2 * stats.t.sf(np.abs(beta / se), dof)
    return beta, se, p


def main():
    T, NT, y = simulate_trios()
    ones = np.ones_like(y)

    # 1) Naive population PGS regression (transmitted only) -> conflates effects.
    b_pop, se_pop, p_pop = ols(np.column_stack([ones, T]), y)
    naive = b_pop[1]

    # 2) Within-family: regress on transmitted AND non-transmitted PGS.
    b_wf, se_wf, p_wf = ols(np.column_stack([ones, T, NT]), y)
    coef_T, coef_NT = b_wf[1], b_wf[2]
    direct = coef_T - coef_NT
    nurture = coef_NT

    print(f"true direct = {BETA_DIRECT:.2f} | true nurture = {BETA_NURTURE:.2f}\n")
    print(f"Naive population PGS effect : {naive:.3f}  (conflates direct + nurture)")
    print(f"Within-family coef(T)       : {coef_T:.3f} ± {se_wf[1]:.3f}")
    print(f"Within-family coef(NT)      : {coef_NT:.3f} ± {se_wf[2]:.3f}  "
          f"(genetic nurture; p = {p_wf[2]:.2e})")
    print(f"  -> recovered DIRECT effect : {direct:.3f}  (true {BETA_DIRECT:.2f})")
    print(f"  -> recovered NURTURE       : {nurture:.3f}  (true {BETA_NURTURE:.2f})")
    print(f"\nNaive overstates direct genetics by {(naive-direct)/direct*100:.0f}%.")

    # ---- figure ----
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    labels = ["Naive PGS\n(population)", "Direct\n(within-family)", "Genetic nurture\n(non-transmitted)"]
    vals = [naive, direct, nurture]
    errs = [se_pop[1], se_wf[1] + se_wf[2], se_wf[2]]
    bars = ax.bar(labels, vals, yerr=errs, capsize=5,
                  color=["#c9a227", "#2a9d8f", "#e76f51"])
    ax.axhline(BETA_DIRECT, ls="--", c="#2a9d8f", lw=1, label=f"true direct {BETA_DIRECT}")
    ax.axhline(BETA_NURTURE, ls=":", c="#e76f51", lw=1, label=f"true nurture {BETA_NURTURE}")
    ax.set_ylabel("standardised effect on outcome")
    ax.set_title("Direct genetic effect vs genetic nurture\n(within-family PGS decomposition)")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width()/2, v + 0.012, f"{v:.2f}", ha="center")
    ax.legend(); fig.tight_layout(); fig.savefig("within_family_results.png", dpi=130)
    print("\nsaved within_family_results.png")

    json.dump({"true_direct": BETA_DIRECT, "true_nurture": BETA_NURTURE,
               "naive_population_pgs": round(float(naive), 4),
               "within_family_coef_T": round(float(coef_T), 4),
               "within_family_coef_NT": round(float(coef_NT), 4),
               "nurture_p": float(p_wf[2]),
               "recovered_direct": round(float(direct), 4),
               "recovered_nurture": round(float(nurture), 4)},
              open("within_family_results.json", "w"), indent=2)
    print("saved within_family_results.json")


if __name__ == "__main__":
    main()
