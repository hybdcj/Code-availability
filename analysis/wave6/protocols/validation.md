# Frozen temporal validation specification (2026-11-14)

English distribution summary of the specification recorded before performance
calculation, after eligibility and event counts had been inspected; not a
prospectively registered original study protocol.
Use the verified four-predictor complete-case risk set for paired model/score
comparison. Main four-year window: 42-54 months, 1,188 participants. Report
2021 (4,041) and 2023 (672) separately. Sensitivities: 45-51 months (1,187),
exclusion of the full development source cohort (1,059), restoration of date
conflicts using official interview months. Freeze all model/score coefficients.
Report AUC and paired AUC difference; continuous-model offset alpha, joint
calibration slope and Brier; fixed integer risk strata. Use 2,000 draws of
2018 communities with replacement and percentile intervals. Keep all members
of sampled communities, weighted by draw multiplicity, without survey weights.
Seed is 20261114 plus sample-ID fingerprint. Plot df=4 cubic B-spline calibration
on the prediction-logit scale over the 1st-99th percentiles; four-year pointwise
bootstrap bands, no forced splitting of tied decile values.
Four-year intercept update expit(lp+alpha) is fitted and evaluated in the same
sample in this stage: apparent performance only, not optimism corrected.
Three/five-year calibration and Brier are cross-horizon diagnostics. Learning
trial and word-list differences, survival, response selection and date proxies
remain limitations. No DCA, IPCW, new Z-score outcome or full coefficient refit.
