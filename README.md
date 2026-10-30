# Mutational Signature Detective

COMP 6721 course project by Anas Farooqui (contact), Hugo Xu and Bao Tran Nguyen.

We estimate how much each known COSMIC mutational signature (for example SBS4 from tobacco or SBS7
from UV light) contributed to a tumor, using the tumor's mutations counted over the 96 standard
substitution contexts. We compare several methods and measure the mutation count below which the
estimates stop being reliable. Our outputs are research estimates, not diagnoses.

## Quick start

```
pip install -r requirements.txt
pytest -q                  # 22 tests
python run_checklist.py    # synthetic demo, about 20 seconds
```

## Real-data run

1. `python gdc_download.py` downloads the open-access TCGA mutation files (LUAD, SKCM, COAD) and
   clinical data into `data/raw/`, checking every file against its GDC checksum.
2. We unpack the GDC reference genome
   ([GRCh38.d1.vd1.fa.tar.gz](https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834))
   into `data/reference/`.
3. `python preprocess_data.py --reference data/reference/GRCh38.d1.vd1.fa` builds the
   96-context mutation counts in `data/processed/`.
4. `python run_checklist.py --signatures data/reference/COSMIC_v3.4_SBS_GRCh38_exome.txt`
   writes `results/run_summary_real.md` and `figures/cosine_vs_mutation_count_real.png`.

The signature list for each cancer type comes from `python derive_signature_subsets.py`: we keep
the signatures that Alexandrov et al. (2020) found in at least 5% of that type's TCGA tumors,
excluding known sequencing artifacts.

## Layout

| Path | Contents |
| --- | --- |
| `sigdet/` | Preprocessing, splits, simulator, methods (B0, M0, M1, M2) and metrics |
| `results/`, `figures/`, `splits/` | Our outputs; `results/week8_progress.md` is the Week 8 write-up |
| `data/processed/`, `data/reference/` | Mutation counts, clinical labels and small reference tables |
| `data/raw/`, genome files | Not committed; rebuilt from `data/manifest.csv` and the steps above |

## Team

- Anas Farooqui: data download and preprocessing, splits, storage.
- Hugo Xu: simulator, classical methods, statistical tests.
- Bao Tran Nguyen: neural network (M3), evaluation and figures.

`TEAM_GUIDE.md` explains every file in the repository and how a run fits together. Limitations and the
plan for Weeks 9 to 13 are in `results/week8_progress.md`.
