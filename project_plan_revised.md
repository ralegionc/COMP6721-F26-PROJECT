# Mutational Signature Detective: revised project plan (draft, 2026-11-10)

Team: Anas Farooqui (40347862, contact), Hugo Xu (40064100), Bao Tran Nguyen (40257379).

This revision replaces the one-page plan (`Mutational Signature Detective Project Plan.pdf`). It keeps
the original goal, methods and checklist, and corrects four things: the owners now match who does each
part of the work, the schedule uses course weeks 8 to 13, the Question 4 answer uses the measured share
of "absent" labels, and the PCAWG cohort is now an optional extension.

## Goal and research questions

We estimate how much each known COSMIC signature (for example SBS4 from tobacco smoke, SBS7 from UV
light) contributed to a tumor, using the tumor's mutation counts over the 96 single-base substitution
contexts. We compare five methods and measure when the estimates can be trusted. Our outputs are
research estimates, not diagnoses.

- **RQ1.** How accurately do non-negative least squares (NNLS), elastic-net NNLS, a maximum-likelihood
  fit and a small neural network recover known exposures?
- **RQ2.** At what mutation count does accuracy break down?
- **RQ3.** How much does accuracy drop when a tumor carries a signature that the method does not know
  about? (Optional extension: does performance hold under a cohort shift?)
- **RQ4.** Is the SBS4 weight higher in smokers with lung adenocarcinoma than in lifelong non-smokers?

## Data, methods and evaluation

**Data.** We use TCGA open-access "Masked Somatic Mutation" MAFs (GDC workflow "Aliquot Ensemble Somatic
Variant Merging and Masking"), filtered to MuTect2-supported calls, on the GRCh38 genome, together with
the GDC clinical records, for lung adenocarcinoma (LUAD), skin melanoma (SKCM) and colon adenocarcinoma
(COAD). The signatures are COSMIC SBS v3.4, exome-normalized. For each cancer type we fit only the
signatures that Alexandrov et al. (2020) found in at least 5% of that type's TCGA exomes, with sequencing
artifacts excluded: 6 for LUAD, 7 for SKCM and 8 for COAD. PCAWG whole genomes are used only if the
optional cohort-shift experiment (E4) goes ahead; in that case both cohorts switch to the MC3 public MAF
so that they share the GRCh37 genome build.

**Ground truth.** True exposures exist only for simulated tumors. We split patients 60/20/20 by their
12-character TCGA patient ID, stratified by cancer type, before calibrating the simulator, and calibrate
it on training patients only.

**Methods.** B0 cohort-mean baseline (every tumor receives the average training exposure), M0 NNLS,
M1 elastic-net NNLS (penalty strength and L1 share tuned on validation data only), M2 multinomial
maximum likelihood fitted with the expectation-maximization (EM) algorithm, with bootstrap intervals,
and M3 a small multilayer perceptron (MLP) trained on simulated tumors. None of M0 to M2 is Bayesian.

**Experiments.** E0 pipeline check; E1 method comparison; E2 accuracy against mutation count (25 to
5,000); E3 robustness to signatures missing from the candidate list; E5 SBS4 against smoking status,
with negative controls. Optional: E4 cohort shift (by tissue source site or PCAWG).

**Metrics.** Exposure cosine similarity, mean absolute error (MAE) on truly present signatures,
per-signature presence F1 and PR-AUC (area under the precision-recall curve), the confusion matrix of
dominant signatures, and the breakdown point N\* (the smallest mutation count with median cosine of at
least 0.90). Methods are compared with paired Wilcoxon signed-rank tests and Holm's correction for
multiple comparisons. Accuracy is reported only to show how much "always absent" would inflate it.

## Schedule and owners

| Week (dates) | Work | Owner | Status |
| --- | --- | --- | --- |
| 8 (Oct 25 to 31) | Download, reference files, signature lists, 96-context catalogs, clinical labels, patient split, storage; first real run with B0, M0, M1, M2 | Anas | Done |
| 9 (Nov 1 to 7) | Simulator and classical methods review; paired Wilcoxon and Holm tests in the pipeline; M1 penalty grid widened | Hugo | Done |
| 9 to 10 (Nov 1 to 14) | M3 MLP trained on simulated training tumors, tuned on validation tumors | Bao Tran | In progress |
| 10 (Nov 8 to 14) | E3 unseen signatures (hold out SBS10a/b and SBS17a); sparse-tumor stress test | Hugo | Planned |
| 10 to 11 (Nov 8 to 21) | E5 extended: more negative controls, presence threshold that grows as mutation count falls | Anas | Planned |
| 11 (Nov 15 to 21) | Full method comparison including M3, statistical tests, M2 bootstrap interval coverage | Hugo | Planned |
| 12 (Nov 22 to 28) | Error analysis (SBS1/SBS5/SBS40a confusion, low-count tumors) and final figures | Bao Tran | Planned |
| 13 (Nov 29 to Dec 5) | Report, reproducibility check from a clean clone, presentation | All | Planned |

Roles: Anas owns data acquisition and preprocessing (everything up to the Week 8 checkpoint) and the
real-data checks; Hugo owns the simulator, classical methods (B0, M0, M1, M2) and statistical tests;
Bao Tran owns M3, evaluation and figures. Everyone reviews pull requests before they are merged.

**Minimum deliverable.** B0, M0, M1, M2 and M3 on all three cancer types with E0, E1, E2 and E5, and the
statistical comparison. E3 is expected; E4 is optional.

**Risks.**
1. The simulator may be too easy, because it generates tumors from the same model the methods assume.
   We report E3 next to E1 and keep the real-data checks of E5.
2. M3 may learn the simulator rather than tumors. We compare it with M0 on the E5 real-data checks too.
3. Absent signatures absorb weight in low-count real tumors (29% of LUAD patients get at least 5% UV
   weight). We test a count-dependent presence threshold in E5.

## Checklist

1. **Is the target clearly defined?** Yes. Each tumor's exposure vector over its cancer type's
   candidate signatures: non-negative shares that sum to 1, known only for simulated tumors.
2. **Can source, patient or object identifiers leak across partitions?** No. We split by patient ID
   before calibrating the simulator, keep all of a patient's samples in one partition, and check the
   split for overlap on every run. If PCAWG is added, donors also present in TCGA are removed.
3. **Is the positive class stated?** Yes. A signature is present at a weight of at least 5%. For the
   smoking test, current or reformed smokers are positive, lifelong non-smokers negative, and patients
   with unknown status are excluded.
4. **Would majority-class accuracy overstate success?** Yes in SKCM and COAD, where 53% and 51% of
   presence labels on the simulated test sets are "absent" (39% in LUAD). The original plan estimated
   about 90%; the final signature lists contain fewer rare signatures. We report F1, PR-AUC and cosine,
   and every method must beat B0.
5. **Can other classifiers reuse precisely the same split?** Yes. The split and the simulated sets are
   saved once with seed 6721 and SHA-256 fingerprints, and every method, including M3, reads the same
   files.
6. **Where will figures and intermediate results be stored?** In the project repository
   (github.com/ralegionc/COMP6721-F26-PROJECT): `results/`, `figures/` and `splits/`. Raw data and
   simulated sets are not committed; they are rebuilt from the download manifest and the seeds.
