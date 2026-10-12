# Pheropsophus-javanus-defense-analysis
> **Table renumbering note (v1.1.0, 2026-10-12):** In the final revision of the manuscript
> (INZ-2026-07-520) the Supplementary Tables were renumbered so that table numbers follow
> their order of first citation. **This repository uses the final, published numbering
> throughout.** If you are comparing against the v1.0.0 archive, the mapping
> (v1.0.0 → final) is: **S24 → S26, S25 → S24, S26 → S25, S33 → S34, S34 → S35,
> S35 → S36, S36 → S33.** All other table numbers are unchanged.

Analysis scripts accompanying the manuscript:

**"Multi-omics evidence for genomic streamlining and co-option of conserved
cellular modules in the explosive chemical defense of the bombardier beetle
Pheropsophus javanus (Dejean, 1825)"**
(resubmission of INZ-2026-07-520, *Integrative Zoology*)

Raw sequencing reads and the genome assembly are deposited at NGDC under
BioProject **PRJCA050879** (GSA: CRA033330; Genome Warehouse: GWHHAIB00000000.1).

## Repository layout

| Folder | Contents | Manuscript results |
|---|---|---|
| `01_gene_family_evolution/` | CAFE5 post-processing and Benjamini-Hochberg corrections | Tables S12-S20 |
| `02_relax_analysis/` | RELAX (HyPhy) batch pipelines, single-copy and multicopy arms | Table S22, Figure S7 |
| `03_transcriptomics/` | DESeq2 analyses and the T. castaneum within-species baseline | Figures 3, Tables S23, S24, S25, S30 |
| `04_genome_surveys/` | GOx genomic survey, OR pseudogenization checks, B. crepitans homolog survey | Tables S7, S33, S36 |
| `05_figures/` | Plotting scripts for the main and supplementary figures | Figures 3A/3B, Figure S7 |

## 01_gene_family_evolution/

- `extract_significant.sh` - extracts nominally significant CAFE5 families and
  their member gene lists/sequences. Usage: `bash extract_significant.sh
  [proteome.fasta]`. Inputs: CAFE5 Base-model outputs (`Base_family_results.txt`,
  `Base_change.tab`, `Orthogroups.txt`) and the P. javanus protein FASTA.
- `CAFE5_BH_correction.py` - Benjamini-Hochberg correction of CAFE5 family- and
  branch-level P values (10,223 tested families). Reproduces the 37 expanded /
  601 contracted families at q < 0.05 (Tables S12-S17). See the script docstring
  for exact input formats.
- `selection_BH_correction.py` - Benjamini-Hochberg correction of the PAML
  selection tests (branch-site PSG arm; branch-model REG arm) across all 2,568
  tested single-copy orthologs, starting from the nominal-significant lists
  (LRT statistics). Reproduces 147 PSGs and 183 REGs at q < 0.05 (Table S20).
  Usage:
  `python3 selection_BH_correction.py --input data/Table_S18_PSGs.tsv --arm PSG --n_tested 2568 --output out.tsv`
  (and likewise for `data/Table_S19_REGs.tsv` with `--arm REG`).
- `data/` - the nominal PSG/REG lists (TSV exports of Supplementary Tables S18
  and S19), the CAFE5 BH q-value table, and the corrected family/gene ID lists
  used by downstream analyses.

## 02_relax_analysis/

- `singlecopy/relax_singlecopy.sh` - RELAX batch pipeline for 42 single-copy
  families (control: stable families). Reuses the species-tree topology, pruned
  per family; helper scripts in `singlecopy/scripts/`.
- `multicopy/relax_multicopy.sh` - RELAX batch pipeline for 241 multicopy
  families, including per-family gene-tree building, codon-alignment hygiene
  (stop-codon sanitization), and run validation; helper scripts in
  `multicopy/scripts/`.
- Family ID lists for each arm sit next to the corresponding script
  (`multicopy/multicopy_241_family_IDs.txt`,
  `singlecopy/controlB_stable_123_singlecopy_IDs.txt`); `data/` holds the
  single-copy ID lists of the other two arms (expanded control, contracted test).
- Required inputs not included (too large / derivable from public data): codon
  alignments per orthogroup and the species tree (see Methods).

## 03_transcriptomics/

- `DESeq2_within_species_gland_vs_headthorax.R` - within-species DESeq2,
  defensive gland vs. head-thorax (n = 3 per tissue). Input: `all_counts.txt`
  (featureCounts matrix; rows = genes, columns = the six samples).
- `DESeq2_cross_species_PJA_vs_TCA.R` - cross-species DESeq2 on reciprocal
  best-hit orthologs. Input: `RBH_with_expression.tsv` (RBH table with
  per-sample expression columns `PJA_expr_1..3`, `TCA_expr_1..2`).
- `TCA_baseline_SRA_pipeline.sh` - download and quantify the T. castaneum
  abdominal gland / non-gland abdomen RNA-seq samples (PRJNA242780;
  SRR1206008, SRR1206010, SRR1206001), single-end.
- `TCA_baseline_validation.R` - difference-of-differences validation of the
  PJA_HEGs against the within-species T. castaneum gland baseline
  (Tables S24, S25).

## 04_genome_surveys/

- `GOx_tblastn_pipeline.sh` - genomic tBLASTn / proteome BLASTp survey for
  glucose-oxidase-like GMC oxidoreductases with 19 insect queries (Table S7).
  Usage: `bash GOx_tblastn_pipeline.sh <genome.fa> [proteins.fa]`.
- `OR_expression_check.R` - normalized expression of the 20 odorant-receptor
  genes from the same `all_counts.txt` matrix (Table S36).
- `OR_CDS_integrity.py` - coding-sequence integrity check (internal stop
  codons / frameshifts). Usage: `python3 OR_CDS_integrity.py cds.fasta`.
- `Bc_presence_survey/` - targeted tBLASTn/BLASTx survey testing whether the
  six principal contracted gene families retain homologs in the Brachinus
  crepitans gland and combined-tissue transcriptomes (Vogel et al., 2025;
  Table S36, Part D). Self-contained; see its own README.

## 05_figures/

- `Figure3A_volcano.R` / `Figure3B_volcano.R` - evidence-tiered volcano plots
  with candidate-gene labels (Figure 3). Inputs in `data/`: gene-annotation and
  baseline-validation tables; the count matrices are as described above.
- `FigS7_RELAX_plot.py` - RELAX three-group comparison figure (Figure S7).
  Usage: `python3 FigS7_RELAX_plot.py data/FigS7_RELAX_input.csv`
  (input identical to Table S22).

## Analyses documented in the Methods (not deposited as scripts)

Standard, parameter-documented pipelines: genome assembly and Hi-C scaffolding,
annotation and BUSCO, OrthoFinder, phylogenomics (IQ-TREE) and molecular dating,
PAML codeml model fitting (its outputs are corrected by
`01_gene_family_evolution/selection_BH_correction.py`), GO/KEGG enrichment
(clusterProfiler), molecular docking (AutoDock4, 50 LGA runs per complex;
AutoDock Vina cross-validation; 2.0 A RMSD pose clustering - parameters in the
Supplementary Methods), and ancestral state reconstruction (phytools
`make.simmap`, ER model selected by AICc, 1,000 simulations; trait coding in
Supplementary Table S4, results in Tables S5-S6).

## Software

R (>= 4.x) with DESeq2, clusterProfiler, EnhancedVolcano, pheatmap, phytools;
Python >= 3.8 (standard library only for the correction scripts; pandas +
matplotlib for plotting); HyPhy (RELAX); BLAST+ >= 2.12; seqkit; CAFE5.




