# A neuronal microexon program controls biomolecular condensation — analysis code

This repository contains the computational analysis code for:

> Lei X., Xia Y., Bi S., Hu Y., Chen Y., Wang Z., Lu Z., Wu X., Tian R. *A neuronal microexon program controls biomolecular condensation* (manuscript).

The code is organized into two pipelines:

- **`microexon_analysis/`** — identification and sequence/phase-separation analysis of conserved neural microexons (Figure 1 and related supplementary figures; also feeds the [NeuronMiniX database](https://neuronminix.med.sustech.edu.cn/)).
- **`imaging_analysis/`** — quality control, feature extraction, and statistical analysis of the paired-isoform high-content imaging screen (Figure 2).

---

## Repository structure

```
microexon_analysis/
├── 1-filter_ASevents.ipynb              # Filter VastDB AS events → conserved neural microexons
├── 2 exon length analysis.ipynb         # Exon length segment classification (3–21 / 24–51 / 54–90 / >90 nt)
├── 3-get_microexon_AAsequence.ipynb     # Exon translation + full-length protein sequences
├── 4-mark miniexon.ipynb                # Mark microexon-encoded residues in protein sequence
├── 5-PICNIC_score.ipynb                 # PICNIC / PSPHunter phase-separation scores
├── 6-get_IDR.ipynb                      # Disorder (IDR) prediction with metapredict
├── 7-get_fuzdrop_score.ipynb            # FuzDrop droplet-promoting probabilities (L and S isoforms)
├── 8-microexon_aa_analysis.ipynb        # Amino acid composition (microexon vs. flanking vs. whole protein)
├── plot-exon-length-distribution.ipynb  # Exon length distribution plot
├── plot-PSI.ipynb                       # PSI heatmap across neural / non-neural samples
├── plot-IDR.ipynb                       # Disorder / IDR plots
├── plot-pDP.ipynb                       # FuzDrop pDP plots (microexon-encoded segment vs. rest)
├── plot-aaComposition.ipynb             # Amino acid composition statistics and plots
├── plot-aaProperties.ipynb              # Amino acid property heatmap
├── plot-phase separation.ipynb          # PICNIC / PSPHunter score distributions (Figure 1G)
└── plot-ASD_PSI.ipynb                   # ASD / Dup15q vs. control PSI plots (Figure 6A)

imaging_analysis/
├── 1. image_QC.ipynb                    # Per-image quality control
├── 2. image_analysis.ipynb              # Per-image feature extraction (81 features)
├── 3. summarize_features.ipynb          # Per-isoform summary + per-gene L-vs-S statistics
└── 4. difference_plot.ipynb             # Volcano plots (Figure 2F–G)
```

---

## microexon_analysis

Numbered notebooks form a sequential pipeline; `plot-*.ipynb` notebooks read the pipeline outputs (primarily `outputs/microexon_final.csv`) and generate figures.

### Pipeline (run in order)

1. **`1-filter_ASevents.ipynb`** — Filters VastDB alternative splicing events for conserved, neural-enriched, in-frame microexons. Criteria: mean PSI > 50 across neural samples (retina excluded), neural/non-neural PSI ratio > 5, cross-species conservation, length divisible by 3, length ≤ 90 nt.
   - Inputs: `source_data/PSI_TABLE-hg38.tab`, `SAMPLE_CONFIG-hg38.tab`, `EVENT_CONSERVATION.tab` (VastDB)
   - Output: `outputs/neuron_up_inframe_miniexon_conserved.csv` (plus `neuron_up_exon_conserved.csv`, `neuron_up_inframe_exon_conserved.csv`)

2. **`2 exon length analysis.ipynb`** — Classifies candidate exons into length segments (3–21, 24–51, 54–90, >90 nt).

3. **`3-get_microexon_AAsequence.ipynb`** — Translates each microexon using CDS phase from GENCODE v47 annotation and retrieves full-length protein sequences from the Ensembl REST API; computes host-protein lengths.
   - Inputs: `source_data/gencode.v47.basic.annotation.gff3`, `EVENT_INFO-hg38.tab`, `PROT_ISOFORMS-hg38.tab` (VastDB)
   - Output: `neuron_up_inframe_miniexon_conserved_aaseq.csv`

4. **`4-mark miniexon.ipynb`** — Marks the microexon-encoded residues with `#` within the full-length protein sequence (with manual curation via `source_data/metadata.csv`).
   - Output: `outputs/miniexon_aa_seq_marked_curated.csv`

5. **`5-PICNIC_score.ipynb`** — Compares phase-separation propensity scores (PICNIC; PSPHunter) of microexon-containing genes against all human genes (KDE plots).
   - Inputs: `source_data/PICNIC-9606-data.csv`, `AllhumanReviewed.txt` (PSPHunter)

6. **`6-get_IDR.ipynb`** — Predicts intrinsic disorder with metapredict; computes exon and protein disorder scores, disordered-domain boundaries, whether the microexon falls within an IDR, and matched random-exon controls.
   - Output: `outputs/miniexon_aa_seq_marked_curated_IDR.csv`

7. **`7-get_fuzdrop_score.ipynb`** — Runs FuzDrop on long (microexon-included) and short (microexon-excluded) isoforms; extracts droplet-promoting probability (pDP) of the exon-encoded segment, its flanks, and the whole protein, plus LLPS probabilities per isoform.
   - Outputs: `outputs/miniexon_aa_seq_marked_curated_IDR_pDP.csv`, per-event profiles in `outputs/fuzdrop/`

8. **`8-microexon_aa_analysis.ipynb`** — Computes amino acid composition and physicochemical descriptors (modlamp) for the microexon-encoded segment, the flanking sequence, and the whole protein.

### Plotting notebooks

| Notebook | Contents |
|---|---|
| `plot-exon-length-distribution.ipynb` | Exon length distribution |
| `plot-PSI.ipynb` | PSI heatmap across neural and non-neural VastDB samples |
| `plot-IDR.ipynb` | Disorder score / IDR localization plots |
| `plot-pDP.ipynb` | FuzDrop pDP of microexon-encoded segment vs. rest of protein |
| `plot-aaComposition.ipynb` | Amino acid composition comparison (chi-square / Fisher exact / Wilcoxon tests) |
| `plot-aaProperties.ipynb` | Amino acid property heatmap |
| `plot-phase separation.ipynb` | PICNIC and PSPHunter score distributions with Mann–Whitney U tests (Figure 1G) |
| `plot-ASD_PSI.ipynb` | PSI of characterized microexons in ASD / Dup15q vs. control postmortem brain (Figure 6A); reads per-event CSVs (`GENE`, `EVENT`, `PSI`, `group`, `LENGTH`) from `source_data/event_csv_BA9_BA41/` (VastDB) |

---

## imaging_analysis

Sequential pipeline analyzing the paired-isoform high-content imaging screen (EGFP-fused long/short isoforms in HEK293T cells, 96-well plates). Expected input layout:

```
images/
└── PROTEIN/
    ├── PROTEIN-L/   # images of the long (microexon-included) isoform
    └── PROTEIN-S/   # images of the short (microexon-excluded) isoform
```

1. **`1. image_QC.ipynb`** — Computes per-image QC metrics (mean foreground intensity, foreground fraction, focus score from normalized Laplacian variance, saturated-pixel fraction) after Triangle-threshold foreground segmentation; applies exclusion thresholds; copies passing images to `images_qc_pass/`.
   - Outputs: `qc_report.csv`, `qc_isoform_summary.csv`, QC summary figures

2. **`2. image_analysis.ipynb`** — Extracts 81 features per QC-passing image (after 2× downsampling and max normalization): intensity statistics (CV, skewness, kurtosis); a granularity spectrum (10 cumulative binary-erosion steps; CellProfiler/Ravkin-style); a Moran's I spatial correlogram (Chebyshev ring lags 1, 3, 5, 10, 20 px) plus a decay term; GLCM texture contrast/correlation (6 distances × 4 angles); and a 50-point foreground-intensity quantile grid.
   - Output: `outputs/per_image_features.csv`

3. **`3. summarize_features.ipynb`** — Summarizes features per isoform (mean/median/SD/MAD) and computes per-gene long-vs-short statistics: Cohen's d (S − L) and two-sided Mann–Whitney U on per-image profile AUCs, plus secondary metrics (cosine distance, Hotelling's T², KS test). Control proteins (EGFP, ATXN1-52Q) are excluded from the gene-level statistics.
   - Outputs: `isoform_summary.csv`, `protein_difference_summary.csv`

4. **`4. difference_plot.ipynb`** — Volcano plots of isoform differences (Cohen's d, S − L; −log₁₀ P, Mann–Whitney U; significance at P < 0.01) for granularity and Moran's I, with the EGFP/ATXN1-52Q control comparison as reference (Figure 2F–G).

---

## Requirements

Analysis was performed with:

- Python 3.12.2
- numpy 2.4.6, scipy 1.17.1, pandas 3.0.3
- scikit-image 0.26.0, scikit-learn 1.9.0
- matplotlib 3.10.9, tifffile 2026.5.15, imagecodecs 2026.5.10
- additionally: seaborn, biopython, requests, metapredict, modlamp, adjustText

External resources: VastDB (splicing events, PSI tables, conservation, protein isoforms), GENCODE v47, Ensembl REST API, FuzDrop, PICNIC, PSPHunter.


## Citation

If you use this code, please cite:

> Lei X., Xia Y., Bi S., Hu Y., Chen Y., Wang Z., Lu Z., Wu X., Tian R. *A neuronal microexon program controls biomolecular condensation* (manuscript).
