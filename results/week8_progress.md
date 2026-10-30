# Week 8 progress checkpoint: Mutational Signature Detective

Team: Anas Farooqui (40347862, contact), Hugo Xu (40064100), Bao Tran Nguyen (40257379).
We ran the pipeline on 2026-10-22. All numbers below come from `results/run_summary_real.md` and the
CSV files next to it, which we produced with
`python run_checklist.py --signatures data/reference/COSMIC_v3.4_SBS_GRCh38_exome.txt`.

## 1. Our data pipeline (done)

| Step | What we did |
| --- | --- |
| Download (`gdc_download.py`) | We downloaded 1,547 open-access TCGA "Masked Somatic Mutation" MAFs (LUAD 618, SKCM 472, COAD 457; 239 MB) and checked each against its GDC MD5 checksum, plus 1,516 clinical records. `data/manifest.csv` lists every file with its SHA-256 checksum. |
| Reference files | We use the GDC GRCh38.d1.vd1 genome and the COSMIC v3.4 SBS GRCh38 exome-normalized signatures (distributed with SigProfilerAssignment, 86 signatures). |
| Preprocessing (`preprocess_data.py`) | We keep single-base substitutions called by MuTect2 in primary or metastatic tumor samples, remove duplicates, and fold them onto the 96 standard contexts. This gives 723,536 mutations in 1,536 tumor samples from 1,451 patients. |
| Genome-build check | None of the 723,536 MAF reference bases disagrees with the GRCh38 genome, so our MAFs and reference match. |
| Clinical labels | In LUAD we have 433 smokers (current or reformed), 75 lifelong non-smokers and 77 patients with unknown status, whom we exclude. Every GDC smoking-status text value maps to a label. |

The median number of mutations per tumor sample is 194 in LUAD, 390 in SKCM and 121 in COAD. Many
tumors have fewer than 100 mutations (LUAD 181 of 613, SKCM 86 of 469, COAD 154 of 454), which is the
range where our breakdown analysis matters.

**Candidate signatures per cancer type.** We derived the signature list for each cancer type from the
per-sample signature attributions for TCGA exomes published with Alexandrov et al. (2020). We keep a
signature for a cancer type when the published fit attributed mutations to it in at least 5% of that
type's tumors and COSMIC does not label it a sequencing artifact (`derive_signature_subsets.py`,
details in `results/signature_subsets.csv`):

- LUAD: SBS1, SBS2, SBS4, SBS5, SBS13, SBS40a
- SKCM: SBS1, SBS5, SBS7a, SBS7b, SBS7c, SBS7d, SBS38
- COAD: SBS1, SBS2, SBS5, SBS13, SBS15, SBS17b, SBS40a, SBS44

The 2020 attributions use COSMIC v3.0 names. COSMIC v3.4 split the old SBS40 into SBS40a, SBS40b and
SBS40c, so we use SBS40a, its closest single successor, in its place.

## 2. Our answers to the six checklist questions (real data)

1. **Target.** Our target is each tumor's exposure vector: the non-negative share of its mutations
   attributed to each candidate signature of its cancer type, summing to 1. We know the true values only
   for simulated tumors.
2. **Identifier leakage.** We split patients 60/20/20 by 12-character TCGA patient ID, stratified by
   cancer type (871 / 290 / 290 patients). The 48 patients with more than one tumor sample keep all their
   samples in one partition. We calibrate the simulator on training patients only.
3. **Positive class.** We call a signature present when its weight is at least 5%. For the smoking
   test, our positive class is current or reformed smokers and our negative class is lifelong
   non-smokers.
4. **Majority-class accuracy.** On our simulated test sets, 39% (LUAD), 53% (SKCM) and 51% (COAD) of
   signature labels are "absent". In SKCM and COAD, predicting "absent" everywhere would score above 50%
   accuracy with zero recall, so we do not report accuracy as a success measure. (Our plan estimated
   about 90% absent labels; our final lists contain fewer rare signatures, so the actual rate is lower.)
5. **Same split for all methods.** We wrote the split (`splits/patient_split_real.json`) and nine
   simulated sets once, with seed 6721 and SHA-256 fingerprints (`splits/sim_fingerprints_real.json`).
   When we reran the pipeline, every result file came out identical.
6. **Storage.** We keep results, figures and splits in `results/`, `figures/` and `splits/` in our
   project repository. We do not commit raw MAFs, the genome or simulated sets; anyone can rebuild them
   from the manifest, the download scripts and the seeds.

## 3. Our baseline and methods (simulated test tumors)

Our simulator draws exposures and mutation counts from distributions fitted to the training patients,
and adds extra-multinomial noise and perturbed signatures. The results below pool mutation counts from
25 to 5,000.

| Cancer | Method | Median cosine | MAE (active) | Macro F1 | PR-AUC | N* |
| --- | --- | --- | --- | --- | --- | --- |
| LUAD | B0 cohort mean | 0.859 | 0.131 | 0.647 | 0.614 | not reached |
| LUAD | M0 NNLS | 0.959 | 0.089 | 0.870 | 0.940 | 50 |
| LUAD | M2 multinomial EM | 0.966 | 0.080 | 0.873 | 0.946 | 50 |
| SKCM | B0 cohort mean | 0.931 | 0.111 | 0.409 | 0.473 | 25 |
| SKCM | M0 NNLS | 0.966 | 0.090 | 0.762 | 0.840 | 50 |
| SKCM | M2 multinomial EM | 0.980 | 0.073 | 0.800 | 0.864 | 25 |
| COAD | B0 cohort mean | 0.913 | 0.100 | 0.444 | 0.492 | 2,500 |
| COAD | M0 NNLS | 0.930 | 0.101 | 0.749 | 0.822 | 250 |
| COAD | M2 multinomial EM | 0.955 | 0.082 | 0.779 | 0.844 | 100 |

N* is the smallest tested mutation count at which the median cosine reaches 0.90. Our M1 (elastic-net
NNLS) gives the same numbers as M0 because validation selected a near-zero penalty (1e-7) for all three
types. Figure: `figures/cosine_vs_mutation_count_real.png` (the M0 curve is hidden under M1).

**Stronger-method progress.** Our M2 is a maximum-likelihood fit of a multinomial model using the EM
(expectation-maximization) algorithm. It beats NNLS on every summary metric in all three types and
lowers the breakdown point for COAD from 250 to 100 mutations. We ran a preliminary paired Wilcoxon
signed-rank test on the same 400 simulated test tumors per type, with Holm correction for three
comparisons. It finds M2's cosine higher than NNLS's in all three types (adjusted p below 1e-6),
although the median gain is small (0.005 to 0.017). We ran this test outside the pipeline and plan
the formal test suite for Week 11.

**Real-data check (E5).** When we apply NNLS to real LUAD tumors, it assigns a median SBS4 (tobacco)
weight of 0.50 to smokers and 0.03 to lifelong non-smokers (424 vs 69 patients with mutation data;
one-sided Mann-Whitney p = 3.7e-21). Among smokers, we find no significant correlation between SBS4
weight and reported pack-years (Spearman rho = 0.08, p = 0.16, n = 348).

## 4. Interpretation and limitations

- **Our simulated accuracy is an upper bound.** We generate simulated tumors from the same signature
  model our methods assume, so a breakdown point of 50 mutations is optimistic. The pilot in our
  project plan found that NNLS needed 500 to 1,000 mutations. We designed the E3 (unseen signature)
  experiment to measure how much accuracy falls when the model is wrong.
- **Exposure cosine flatters the baseline when one signature dominates.** B0, which predicts the
  training-cohort average for every tumor, reaches a median cosine of 0.93 in SKCM because SBS7a
  dominates most melanomas. Its macro F1 of 0.41 shows that it does not identify which signatures are
  present, so we cannot rely on cosine alone.
- **Signature leakage on real data.** In a negative control, we added the UV signatures SBS7a to SBS7d
  to the LUAD fit. They receive at least 5% weight in 29% of LUAD patients, although we do not expect
  lung tumors to carry UV damage. The reverse control (SBS4 added to SKCM) gives 9%. Signatures that
  should be absent absorb part of each tumor's mutations, mostly in low-count tumors, so our 5%
  presence threshold is too low for real tumors with few mutations.
- **Data limitations.** Open-access MAFs remove some low-frequency variants, which lowers counts in
  low-burden tumors. Our candidate lists exclude signatures seen in under 5% of tumors (for example the
  POLE signatures SBS10a and SBS10b in about 2% of colorectal tumors), so we will fit tumors driven by
  them poorly. The published attributions predate COSMIC v3.4, so newer signatures (for example SBS92
  and SBS88) are not candidates.

## 5. Our plan for Weeks 9 to 13 and risks

| Week | Tasks |
| --- | --- |
| 9 | Train the M3 neural network on simulated tumors; widen M1's penalty grid and tune M1 and M2 on validation data |
| 10 | E3 robustness to unseen signatures (for example holding out SBS10a/b and SBS17a); E5 with more negative controls and a count-dependent presence threshold |
| 11 | Paired Wilcoxon tests with Holm correction inside the pipeline; optional E4 cohort shift by tissue source site |
| 12 | Error analysis: confusion between similar signatures (SBS1/SBS5/SBS40a), low-count tumors; final figures |
| 13 | Report, reproducibility check from a clean clone, presentation |

Risks: (1) Our simulator may be too easy, which inflates every method's accuracy; we will add
model-mismatch noise and report E3 alongside E1. (2) M3 may only learn the simulator; we will compare
it with M0 on the E5 real-data checks, not only on simulations. (3) The 3 GB genome takes disk space;
after preprocessing we only need the mutation counts in `data/processed/`.

## References

Alexandrov, L. B., et al. (2020). The repertoire of mutational signatures in human cancer. *Nature*, 578, 94-101.
COSMIC Mutational Signatures v3.4, https://cancer.sanger.ac.uk/signatures/
