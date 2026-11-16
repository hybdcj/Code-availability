# Intercept update internal-validation specification (2026-11-16)

English distribution summary of the specification recorded before this analysis.
Use saved four-year predictions: 1,188 participants, 269 events, 212 communities.
Freeze all original coefficients. With seed 20261116, draw 212 communities with
replacement 10,000 times. Retain every participant in each sampled community
and weight by community multiplicity. Sample sizes vary; Brier is averaged
over participants, not equally over communities.
Refit only the offset intercept in each draw. Evaluate that updated model on
the bootstrap training sample and the entire original sample. Since lower
Brier is better, correction = mean(original-test loss - bootstrap-training
loss); corrected Brier = full-sample apparent Brier + correction. Compare with
the original frozen model's validation Brier. Do not silently replace failures.
Save replicate losses, counts and provenance. Report point estimates and the
Monte Carlo SE of the correction, not a sampling CI. No nested bootstrap or
cross-validation. Original-sample testing is neither independent nor out-of-bag.
Verify the full fit and first 50 draws independently with R glm, loss arithmetic,
and unchanged AUC. This validates the fixed intercept-updating procedure, not
earlier predictor/outcome/method selection, and is not external validation.
