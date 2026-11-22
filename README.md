# Solubility Audit

Can six molecular descriptors predict how well an organic compound dissolves in water? And what do the errors teach us about what the models miss?

This project connects a Chem 81 recrystallization question to a public structure–property dataset. It benchmarks two small models and audits their errors using hydrogen-bond donor counts and a phenol–anisole comparison.

## Measured benchmark

**1,128 ESOL records · 5 held-out folds · 6 input descriptors · 2 fitted models**

| Model | MAE (logS) | RMSE (logS) |
|---|---:|---:|
| Training-mean baseline | 1.660 | 2.097 |
| Ridge regression | 0.929 | 1.189 |
| 8-nearest neighbors | **0.709** | **0.938** |

The nearest-neighbor model lowered MAE by **57.3%** relative to the training-mean baseline. Each data row was scored once while held out. Identical source SMILES strings stayed in the same fold. Other representations of the same molecule and shared scaffolds can occur across folds, so these results estimate performance within ESOL under this split.

### Where the model missed

| Hydrogen-bond donor count | Records | Ridge MAE | 8-nearest-neighbor MAE |
|---|---:|---:|---:|
| 0 | 635 | 1.012 | 0.770 |
| 1 | 307 | 0.751 | 0.592 |
| 2+ | 186 | 0.938 | 0.696 |

In ESOL, phenol measured **0.00 logS** and anisole measured **−1.85 logS**. The nearest-neighbor model predicted **−0.67** and **−2.31**, a gap of **1.64** compared with **1.85** measured. Phenol's O–H group can donate a hydrogen bond and anisole's methoxy group lacks that donor. Methylation also changes molecular weight and other properties, so this two-compound comparison cannot isolate a causal effect. Both compounds share a benzene scaffold and were scored in different folds. The separate scaffold experiment groups them together when evaluating the pipeline.

## Reproduce

Requires Node.js 22 for the verified benchmark. Run from the project root.

```sh
node benchmark.mjs
```

The command downloads the public Delaney CSV and writes summary.json and predictions.csv under results/. To reuse a local copy, run `node benchmark.mjs esol.csv`. See [the verified summary](benchmark-results.json) for full precision, fold sizes, input hash and the five largest ridge errors.

For the RDKit experiment, Python 3.11 and the pinned RDKit version regenerate descriptors from SMILES and add Bemis–Murcko scaffold groups. All ring-free compounds share one group. The large groups may produce uneven test folds. The [RDKit scaffold benchmark](https://github.com/nisoo0411/solubility-audit/actions/runs/44654418835) evaluated 1,128 records across 269 scaffold groups using five-fold cross-validation. Each test fold held out entire scaffold groups.

| Model | MAE (logS) | RMSE (logS) |
|---|---:|---:|
| Training-mean baseline | 1.762 | 2.210 |
| Ridge regression | 1.009 | 1.293 |
| 8-nearest neighbors | **0.988** | **1.232** |

The nearest-neighbor model had 44% lower MAE than the scaffold benchmark baseline. Its MAE increased from 0.709 in the exact-SMILES benchmark to 0.988 in the scaffold benchmark. Both the grouping method and the descriptor source changed between runs.

```sh
python -m pip install -r requirements.txt
python rdkit_descriptors.py
node benchmark.mjs rdkit-esol.csv results/scaffold
```

The scaffold experiment changes both the split and descriptor source. Compare the results as two evaluation protocols. To isolate split effects, use the recomputed CSV in both runs and remove only its scaffold column for the exact-SMILES split.

## Methods and data

- Target is measured base-10 log solubility in mol/L. A one-unit prediction error corresponds to a tenfold concentration ratio.
- Inputs are minimum atom degree, molecular weight, donor count, ring count, rotatable-bond count and topological polar surface area. The verified run used these precomputed columns supplied with ESOL. The supplied ESOL prediction column is excluded.
- A deterministic group assignment places each entire group into the smallest fold, with hashed group names breaking size ties. Each model standardizes descriptors from its training fold only. Ridge uses a 0.01 penalty on standardized coefficients in the mean squared loss. Nearest neighbors use Euclidean distance and the mean of 8 neighbors. Settings are fixed for the benchmark and no hyperparameter search is reported.
- Solubility can depend on pH, temperature and solid form. This benchmark does not use those variables and the donor slices are descriptive comparisons.
- The reported figures come from a run against the public measurement data.

Dataset: [DeepChem Delaney ESOL CSV](https://github.com/deepchem/deepchem/blob/master/datasets/delaney-processed.csv). Original paper: J. S. Delaney, *ESOL: Estimating Aqueous Solubility Directly from Molecular Structure* (2004), [DOI 10.1021/ci034243x](https://doi.org/10.1021/ci034243x).
