# Rice immune gene prioritization with PU learning

Prioritizing candidate immune genes in rice (chitin-triggered immunity) with positive-unlabeled learning on sequence-derived features, followed by co-expression network validation.

**My contribution:** I was the sole author. I designed and implemented the full analysis: RNA-seq processing and DESeq2, feature engineering, PU bagging, evaluation, network analysis and figures (M.Tech thesis, Pondicherry University, supervised by Dr. P.T.V. Lakshmi).

## Method in brief
1. **Positives:** 394 DEGs (DESeq2, PRJNA1062417). Everything else is *unlabeled*, not negative.
2. **Features per gene:** promoter GC, W-box and GCC-box motif counts, CDS GC and length, protein length, cysteine count and density, signal-peptide probability.
3. **PU bagging:** 10 rounds; each draws as many unlabeled genes as positives as temporary negatives, trains a Random Forest (300 trees, seeds 0-9) and scores all genes. `PU_score` is the mean probability; `PU_std` its spread.
4. **Checks:** OOB accuracy, stability across rounds, feature importance, Mann-Whitney tests, mutual information.
5. **Network step:** co-expression network (Pearson, 233 nodes / 8,177 edges); PU score is combined with degree centrality to give four non-DEG candidates.

## Layout
- `notebooks/rice_immune_PU_learning.ipynb`: the analysis, top to bottom
- `results/tables/`: top integrated candidates, feature importances
- `results/figures/`: PCA, volcano, MA, heatmap, PU score distribution, candidate subnetwork
- `DATA.md`: input files and where to get them
- `requirements.txt`

## Reproducibility notes
- Run order is top to bottom after placing the inputs from `DATA.md`.
- Forest seeds are fixed. The exploratory sampling in Section 4 uses unseeded `np.random.choice`, so those numbers vary slightly; the PU bagging in Section 6 is fully seeded.
- Limitations: positives are DEGs from one experiment, so the model learns what resembles DEGs by sequence features; candidates are hypotheses and have no wet-lab validation.

## Full workflow (original project)
RNA-seq (PRJNA1062417) -> FastQC / MultiQC -> Trimmomatic -> HISAT2 -> samtools -> featureCounts (Subread) -> DESeq2 (R) -> feature engineering -> PU learning -> gene prioritization -> co-expression network validation.

The upstream steps (QC to DESeq2) were run with the tools above; the DEG table they produced (`Rice_deg_with_ID.csv`) is the input to the notebook. The notebook covers feature engineering through prioritization.

## Key results
- 394 immune-responsive DEGs
- Genome-wide feature matrix for about 35,000 unlabeled genes
- PU-learning model trained on the immune genes, ranked by stability across 10 rounds
- Network validation: 233 nodes, 8,177 edges
- 4 non-DEG immune candidates prioritized (PU score plus network centrality)

## Software
Python 3.10 (scikit-learn, pandas, numpy, biopython, networkx, seaborn, matplotlib), R with DESeq2, plus the command-line tools listed above.
