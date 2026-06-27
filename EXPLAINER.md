# What this demo does, in plain language

One small program (`within_family_demo.py`) answering the question at the heart of
social-science genetics. Runs in seconds (`numpy`, `scipy`, `matplotlib`).
Data is **simulated** — the point is to show the method is correct, then swap in
real trio data.

## The question (plain)

A "polygenic score" (PGS) adds up many DNA variants to predict an outcome like
educational attainment. It predicts well — but *why*? Two reasons get mixed up:
1. **Direct (nature):** the child's own genes affect the child.
2. **Genetic nurture (nurture):** the *parents'* genes shaped the home
   environment, which affects the child.

A normal PGS analysis cannot tell these apart and so overstates "nature".

## What it does

For each child it splits the parents' DNA into:
- **Transmitted** alleles = the child's own genotype.
- **Non-transmitted** alleles = the half the parents carry but did *not* pass on.

Non-transmitted alleles are **not in the child**, so they cannot act directly —
yet if they still predict the child's outcome, that can only be through the
environment the parents created. Which allele is transmitted is decided at random
by nature (Mendelian segregation) — a built-in natural experiment. Then a simple
regression of the outcome on the transmitted (T) and non-transmitted (NT) scores:

```
coef(NT) = genetic nurture
coef(T)  = direct + nurture
direct   = coef(T) − coef(NT)
```

## Result (simulated, reproducible)

| Estimate | Value | Truth |
|---|---|---|
| Naive population PGS | **0.45** | (mixes both) |
| Direct effect (within-family) | **0.31** | 0.30 |
| Genetic nurture (non-transmitted) | **0.13**, p ≈ 7e−16 | 0.15 |

The naive PGS **overstates the direct genetic effect by ~42%**. The within-family
design recovers the true direct effect and proves genetic nurture is real.

## Why it matters

This transmitted/non-transmitted design (Kong et al., Science 2018) is the core
method of the group's research. It is also an **instrumental-variables /
natural-experiment** identification — exactly the econometric toolkit — applied to
genetics.

## ML / statistics angle

Polygenic scores are a linear genomic predictor; the decomposition is OLS with
an exogenous regressor (NT) created by random Mendelian transmission. The
identification (not the prediction accuracy) is the contribution — a causal-
inference framing rather than a black-box model.

## Honest limitations

- **Simulated data.** The numbers prove the method and code are correct, not a
  real finding. Real use: compute T and NT polygenic scores from phased trio
  genotypes + external GWAS weights (e.g. educational-attainment GWAS, Lee 2018),
  then run the same regression.
- Single-population, no assortative mating, no population stratification modelled
  — kept deliberately simple to show the core identification cleanly.
