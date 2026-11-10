# Team guide

This guide explains how to get the project running, how the pieces fit together, and what every file in
the repository is for. `results/week8_progress.md` describes where we stand scientifically.

## 1. What the project does, in one paragraph

Every tumor carries mutations, and different causes leave different patterns. COSMIC catalogues these
patterns as "signatures" (for example SBS4 from tobacco smoke, SBS7a from UV light). Each signature is a
probability profile over the 96 possible single-base mutation types, where a type is the change itself
(such as C>T) plus the base on each side (for example `A[C>T]G`). We count a tumor's mutations in those
96 bins and estimate what share came from each signature: the tumor's **exposure vector**. Real tumors
have no known answer, so we test our methods on **simulated tumors** whose true exposures we set
ourselves, then run sanity checks on real tumors (for example, smokers' lung tumors should show more
SBS4). We also find the **breakdown point N\***, the smallest mutation count at which the estimates are
still good (median cosine similarity of at least 0.90 between true and estimated exposures).

## 2. Getting started (about 5 minutes)

```
git clone https://github.com/ralegionc/COMP6721-F26-PROJECT.git
cd COMP6721-F26-PROJECT
pip install -r requirements.txt
pytest -q
python run_checklist.py --signatures data/reference/COSMIC_v3.4_SBS_GRCh38_exome.txt
```

- `pytest -q` should report 28 passed.
- The last command reruns the full real-data experiment in under a minute. It does **not** need the raw
  TCGA downloads or the 3 GB reference genome: their processed output is already in `data/processed/`.
- Afterwards, `git status` should show no changed files under `results/`, which confirms that the run
  reproduced our numbers exactly.

**Always run `run_checklist.py` with a flag:**

| Command | What it does |
| --- | --- |
| `python run_checklist.py --signatures data/reference/COSMIC_v3.4_SBS_GRCh38_exome.txt` | Real tumors, real COSMIC signatures. Writes the `*_real` outputs. |
| `python run_checklist.py --demo` | Fully synthetic cohort and fake signatures, for quick tests. Writes the `*_demo` outputs. |
| `python run_checklist.py` (no flag) | **Avoid.** It loads the real tumors but uses the fake demo signatures, and overwrites the real results with meaningless numbers. |

## 3. How a run flows

```
 GDC website ──gdc_download.py──> data/raw/  (not in Git; listed in data/manifest.csv)
                                     │
             reference genome ──preprocess_data.py──> data/processed/*_sbs96.csv, clinical.csv
                                                         │
 COSMIC signatures (data/reference/) ────────────────────┤
                                                         ▼
                                                  run_checklist.py
   1. split patients 60/20/20 (once)          ──> splits/patient_split_real.json
   2. fit NNLS to the real training tumors and calibrate the simulator on them
   3. simulate train/val/test tumors (once)   ──> data/sims/ (not in Git), splits/sim_fingerprints_real.json
   4. fit B0, M0, M1, M2 on train/val; score on the simulated test tumors
   5. paired Wilcoxon tests between methods, Holm-corrected
   6. smoking test and negative controls on real LUAD/SKCM tumors
   7. write tables, figure and summary        ──> results/, figures/
```

Steps 1 and 3 happen only the first time. Later runs reload the saved files and check their SHA-256
fingerprints, so every method is always scored on exactly the same split and the same test tumors.

## 4. Every file and what it does

### Top-level scripts (the entry points)

| File | What it does | When to run it |
| --- | --- | --- |
| `run_checklist.py` | The main experiment. Loads the data, makes or loads the split, simulates tumors, fits every method, computes all metrics, runs the smoking test, and writes `results/` and `figures/`. Also answers the instructor's six checklist questions in the run summary. | Yes, after any change to methods or metrics. |
| `run_method_tests.py` | Reruns only the paired Wilcoxon comparisons from the saved per-tumor scores (`--demo` for the demo run) and rewrites `results/method_tests_<source>.csv`. | Only when changing the tests; `run_checklist.py` already runs them. |
| `gdc_download.py` | Downloads the 1,547 open-access TCGA mutation files (MAFs) and clinical records from the GDC API, checks each file's MD5 checksum, and writes `data/manifest.csv`. `--dry-run` only lists what it would download. | Not needed; its processed output is in `data/processed/`. |
| `preprocess_data.py` | Turns the raw MAFs into mutation counts per tumor over the 96 contexts, using the reference genome to look up neighbouring bases, and builds the clinical table. Needs `--reference <genome .fa>`. | No, unless we change a preprocessing filter. |
| `derive_signature_subsets.py` | Decides which signatures we fit for each cancer type, using published TCGA results (kept if seen in at least 5% of that type's tumors, artifacts excluded). Prints the lists that go into `sigdet/config.py`. | Only if we change the 5% rule. |

### `sigdet/`: the library the scripts use

| File | What it contains |
| --- | --- |
| `config.py` | All fixed settings: random seed (6721), split fractions, 5% presence threshold, 0.90 breakdown threshold, the mutation counts tested (25 to 5,000), the 96 context names, artifact signatures, and the signature list per cancer type. **Changing anything here changes every result.** |
| `methods.py` | The estimation methods, all with the same `fit(train, val)` / `predict(counts)` interface: **B0** (always predicts the training average; our "dumb" baseline), **M0** NNLS (non-negative least squares), **M1** elastic-net NNLS (penalty strength and L1 share chosen on validation data from a 9 x 3 grid, by MAE on present signatures), **M2** multinomial maximum likelihood via the EM algorithm, with bootstrap intervals. **M3 goes here.** |
| `metrics.py` | How we score methods: cosine similarity per tumor, mean absolute error on truly present signatures, per-signature precision/recall/F1/PR-AUC, the confusion matrix of dominant signatures, the breakdown point N\*, the "always absent" check that shows why accuracy is misleading, and the paired Wilcoxon tests with Holm correction that compare methods (`compare_methods`). |
| `simulate.py` | The simulator. `calibrate()` learns typical exposures and mutation counts from **training patients only**; `simulate()` draws tumors with known exposures, adding extra noise and slightly distorted signatures so the test is not too easy. `n_active=(1, 3)` instead gives sparse tumors with 1 to 3 random active signatures. |
| `splits.py` | The patient-level split. Extracts the 12-character patient ID from TCGA barcodes, splits patients (not samples) stratified by cancer type, checks that no patient appears in two partitions, and saves or loads the split with a tamper-check fingerprint. |
| `target.py` | Defines what we predict (the exposure vector) and validates it (non-negative, sums to 1). Also loads the COSMIC signature file and picks each cancer type's subset. |
| `labels.py` | Defines the positive classes: a signature is "present" at 5% or more; for smoking, current or former smokers are positive, lifelong non-smokers negative, unknown excluded. Maps the GDC's smoking-status text to these labels. |
| `realdata.py` | Checks on real tumors, which have no ground truth: the smoking test (is SBS4 higher in smokers?) and negative controls (signatures that should be absent, such as UV in lung tumors). |
| `clinical.py` | Reads the GDC clinical records and extracts smoking status and pack-years per patient. |
| `preprocess.py` | The MAF-to-counts pipeline used by `preprocess_data.py`: keeps single-base substitutions called by MuTect2, removes duplicates, finds each mutation's context, and stops if the genome build is wrong. |
| `reference.py` | A small reader for the 3 GB genome file that jumps straight to a position instead of loading everything. |
| `storage.py` | Folder layout, checksums, the download manifest, and saving/loading simulated sets and result tables. |
| `__init__.py` | Marks `sigdet` as a package and maps each checklist question to its module. |

### `tests/`

| File | What it checks |
| --- | --- |
| `test_checklist.py` | One or more tests per checklist question: valid exposure vectors, no artifact signatures, patient-level split with leak detection, simulator uses training patients only, positive-class definitions, accuracy inflation, NNLS beating B0, breakdown point, split reproducibility and tamper detection, storage layout. |
| `test_preprocess.py` | Preprocessing on a tiny fake genome: context counting and strand folding, wrong-genome detection, caller-column handling, clinical parsing, and that M1/M2 recover known exposures. |
| `test_methods_and_stats.py` | Holm correction against a hand computation, Wilcoxon against SciPy, tumor-by-tumor pairing in `compare_methods`, per-tumor MAE, M1's full grid search, and the sparse simulator option. |

We run them with `pytest -q`, and new code comes with new tests.

### `data/`

| Path | Contents |
| --- | --- |
| `data/processed/LUAD_sbs96.csv`, `SKCM_sbs96.csv`, `COAD_sbs96.csv` | **The main input.** One row per tumor sample (TCGA barcode), 96 columns of mutation counts. 613, 469 and 454 samples. |
| `data/processed/clinical.csv` | One row per patient: cancer project, raw smoking status, pack-years, and our smoker label (1, 0 or empty). |
| `data/processed/qc_*.json` | How many mutations survived each preprocessing filter, per cancer type. |
| `data/reference/COSMIC_v3.4_SBS_GRCh38_exome.txt` | The 86 COSMIC v3.4 signatures (96 rows each), normalized for exome data. Passed to `--signatures`. |
| `data/reference/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv` | Published signature attributions for TCGA tumors (Alexandrov et al. 2020); input to `derive_signature_subsets.py`. |
| `data/manifest.csv` | Every raw downloaded file with its size and SHA-256 checksum, so the raw data can be re-downloaded and verified. |
| `data/raw/`, `data/sims/`, genome files | **Not in Git.** Raw downloads, simulated tumors (rebuilt automatically from the seed on the first run) and the 3 GB genome. |

### `splits/`

| File | Contents |
| --- | --- |
| `patient_split_real.json` | Which patients are in train (871), validation (290) and test (290), with the seed and a fingerprint. **Never edit or regenerate it.** |
| `sim_fingerprints_real.json` | Checksums of the nine simulated sets (3 cancer types x train/val/test). |
| `*_demo.json` | The same for the synthetic demo. |

### `results/` (all written by `run_checklist.py`, except the `.md` write-ups listed first)

Every file exists twice: `_real` (real data, the one that matters) and `_demo` (synthetic).

| File | Contents |
| --- | --- |
| `preliminary_report_draft.md` | Draft of the final report (abstract to references) for lab discussion. **Start here.** |
| `progress_report.md` | Status report: work done so far and Weeks 10 to 13 plan with owners. |
| `week8_progress.md` | Our Week 8 write-up: data pipeline, answers to the six questions, results, limitations and plan. |
| `run_summary_real.md` | Human-readable summary of the latest run: checklist answers, results table, breakdown points, smoking test. |
| `summary_real.csv` | One row per cancer type and method: median cosine, MAE, macro/micro F1, PR-AUC, accuracy. |
| `per_tumor_cosine_real.csv` | Cosine and MAE (present signatures) for every simulated test tumor, with its mutation count and a `tumor` position shared across methods. Input for paired tests and plots. |
| `method_tests_real.csv` | Paired Wilcoxon test for every pair of methods, per cancer type, pooled (`n` = all) and at each mutation count, with Holm-adjusted p-values. |
| `breakdown_point_real.csv` | N\* per method and cancer type. |
| `presence_per_signature_real.csv` | Precision, recall, F1 and PR-AUC for each signature, method and cancer type. |
| `confusion_<TYPE>_<method>_real.csv` | For M0 and M2: how often each true dominant signature was estimated as each signature. |
| `majority_class_check_real.csv` | What an "always absent" guess would score, per cancer type. |
| `e5_smoking_real.csv` | Smoking test: SBS4 in smokers vs never-smokers, Mann-Whitney test, pack-years correlation. |
| `e5_negative_controls_real.csv` | Weight given to signatures that should be absent. |
| `signature_subsets.csv` | Why each signature was or was not selected per cancer type. |

### `figures/`

| File | Contents |
| --- | --- |
| `cosine_vs_mutation_count_real.png` | Median cosine against mutation count for each method and cancer type, with the 0.90 line that defines N\*. |
| `e5_sbs4_smoking_real.png` | Estimated SBS4 weight per LUAD patient, lifelong non-smokers against smokers. |

### Other files

| File | Contents |
| --- | --- |
| `README.md` | Short project overview and run instructions. |
| `requirements.txt` | Python packages to install. |
| `.gitignore` | Keeps raw data, simulations and the genome out of Git. |
| `COMP 6721 Course Project Proposal.pdf`, `Mutational Signature Detective Project Plan.pdf` | Our original proposal and one-page plan. |

## 5. Where new work goes

**Classical methods and statistics**
- M1's search grid is `M1_ALPHAS` and `M1_L1_RATIOS` at the top of the M1 class in `sigdet/methods.py`.
- New statistical tests go in `sigdet/metrics.py` and are called from `run_checklist.py`, so they always
  run on the same saved test tumors. `compare_methods` pairs methods through the `tumor` column.
- Simulator variants (for example sparse tumors via `n_active`) are options of `simulate()` in
  `sigdet/simulate.py`, not a separate simulator.

**Neural network (M3), evaluation and figures**
- M3 is a class in `sigdet/methods.py` with a `name` attribute, `fit(train, val=None)` and
  `predict(counts)`. `train` and `val` are dicts with `"counts"` (tumors x 96) and `"exposures"`
  (tumors x signatures). `predict` must return non-negative rows that sum to 1.
- Adding it to the `methods = [...]` list in `run_checklist.py` scores it on exactly the same test
  tumors as the other methods, and it then appears in every results file and the figure.

**Ground rules**
- Nothing in `splits/` is edited or regenerated, and `sigdet/config.py` changes only by team agreement.
- `data/sims/` and `data/raw/` stay out of Git (`.gitignore` already blocks them).
- Committed results are the `_real` files from a run with `--signatures`. Outputs go to `results/`
  and `figures/`; `data/processed/` holds only the real-data inputs.
- Scripts build paths from `sigdet.storage.ensure_dirs()`, so they work from any folder.

## 6. How we use Git

1. Create a branch per topic: `git checkout -b name/short-topic` (for example `m3-mlp`).
2. Commit, run `pytest -q`, then `git push -u origin <branch>`.
3. Open a pull request on GitHub for review before merging into `main`.
4. Before starting new work: `git checkout main` and `git pull`.
