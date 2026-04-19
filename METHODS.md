# Materials and Methods

## Study design

We performed a retrospective pan-cancer analysis to evaluate whether concordant missense pathogenicity predictions identify clinically informative subsets of HRR VUS carriers, using an integrative framework complemented by LOH and genomic scar characterization. AlphaHRD combines computational pathogenicity prediction (Layer 1) with allele-specific loss of heterozygosity (Layer 2) and genomic instability scores (Layer 3). To test whether predictor concordance sharpens the clinical signal beyond any single tool, we benchmarked AlphaMissense against three additional missense pathogenicity predictors — REVEL, CADD, and PrimateAI — within the same cohort. All survival analyses were performed at the patient level; when a patient carried multiple HRR missense variants, classification was based on the maximum pathogenicity score across variants for each predictor. All data were obtained from public repositories; institutional review board approval was not required. We followed the REMARK guidelines for tumor marker studies [39] and TRIPOD+AI for prediction model reporting [40]. A pre-submission code audit and complete provenance trail are documented in the repository (AUDIT.md).

## Data sources

Somatic mutations were obtained from 31 TCGA PanCancer Atlas studies via the cBioPortal REST API (accessed February 16, 2026) [26,27]. These derive from the MC3 consensus call set, which applies seven mutation callers to whole-exome sequencing across 10,295 tumors [29]. Clinical endpoints included overall survival with event coding verified per study.

AlphaMissense scores were retrieved from Zenodo (doi:10.5281/zenodo.8208688) as per-protein TSV files [14]. Scores from REVEL [REVEL ref], CADD (phred-scaled) [CADD ref], and PrimateAI [PrimateAI ref] were obtained through the myvariant.info batch API (GRCh37 assembly) [myvariant ref]. ClinVar annotations (variant_summary.txt.gz) were downloaded from the NCBI FTP server on February 16, 2026.

Allele-specific copy number segments were retrieved from the Genomic Data Commons using the ASCAT3 workflow [18], providing major and minor allele copy numbers per segment. Pre-computed HRD scores (LOH count, TAI, and LST) came from the GDC PanCan-DDR-2018 supplementary archive, originally reported by Knijnenburg et al. for 9,125 TCGA samples [20]. We verified that HRDsum = LOH + TAI + LST for all records (maximum discrepancy: zero). The threshold of HRDsum ≥ 42 follows Telli et al. [19]. Mutational signatures used the MC3 MAF with SigProfilerAssignment v0.1.6 and COSMIC v3.4 references [32,21]. BRCA1 promoter methylation (probes cg13601799, cg08047457, cg19531713) was accessed via UCSC Xena [30]; immunogenomic features from Thorsson et al. supplementary Table S1 [31].

## Gene panel

Twenty-five HRR and DNA damage response genes were selected in three tiers based on clinical evidence for PARP inhibitor sensitivity: PROfound Cohort A (BRCA1, BRCA2, ATM), PROfound Cohort B (PALB2, BRIP1, BARD1, CDK12, CHEK1, CHEK2, FANCL, RAD51B–D, RAD54L), and an expanded Fanconi anemia/DDR panel (FANCA, FANCC, FANCD2, FANCE–G, NBN, MRE11, RAD50, ATR, ATRX) [4,41,42]. Only missense variants were included, as all four predictors are designed for this variant class.

## Layer 1: Pathogenicity classification and benchmarking

Each variant was mapped to its UniProt accession and annotated with AlphaMissense scores (554,363 pre-computed values across 25 genes; match rate 99.4%). Published categorical thresholds were applied without optimization: AlphaMissense likely pathogenic > 0.564 [14], REVEL ≥ 0.5 [REVEL ref], CADD ≥ 25 phred [CADD ref], PrimateAI ≥ 0.803 [PrimateAI ref]. Global score coverage across all 4,301 variants was 99.4% for AlphaMissense, 87.7% for CADD, 87.6% for REVEL, and 86.8% for PrimateAI. Within the germline-enriched survival subset (n = 911), all four predictors had comparable patient-level coverage. ClinVar concordance was assessed using Cohen's κ with 2,000-iteration bootstrap 95% confidence intervals (seed = 42).

## Layer 2: Loss of heterozygosity

For patients with at least one AlphaMissense-pathogenic variant, ASCAT3 allele-specific segments were intersected with gene coordinates (GRCh38). LOH was defined as minor copy number equal to zero at the variant locus, capturing both hemizygous deletion and copy-neutral LOH. Patients were classified as biallelic (variant plus LOH) or monoallelic (variant present, wild-type allele retained). Coverage reached 632 of 716 patients (88.3%).

## Layer 3: Genomic scar characterization

HRDsum was compared between biallelic and monoallelic groups by Mann–Whitney U test; the proportion exceeding ≥ 42 was compared by Fisher's exact test. SBS3 mutational signature exposures provided an independent measure of HRR deficiency. Patients were assigned to five biological characterization tiers — TRUE_HRD (biallelic + HRDsum ≥ 42), PROBABLE_HRD (biallelic + HRDsum 33–41), BIALLELIC_NO_SCAR (biallelic + HRDsum < 33), MONOALLELIC, and BENIGN — used for genomic characterization rather than as the primary prognostic classifier.

## Survival analysis

Overall survival was the primary endpoint. The primary model was a stratified Cox regression using tumor type as the stratification variable, with AlphaMissense pathogenicity status as the predictor. Per-tumor hazard ratios were synthesized using fixed-effect and REML random-effects meta-analysis with the Hartung–Knapp adjustment [35,36].

The key secondary analysis tested whether pathogenicity reclassification — and cross-predictor concordance specifically — adds prognostic value within germline-enriched genes. We restricted to patients carrying variants in BRCA1, BRCA2, ATM, PALB2, and CHEK2, where the majority of pathogenic variants are germline in origin (Huang et al. [23]). Each predictor was tested individually (pathogenic vs. benign within these five genes). For the concordance analysis, patients were classified into four groups based on AlphaMissense and REVEL calls jointly: both pathogenic, AM-only pathogenic, REVEL-only pathogenic, and both benign (reference). A gene-stratified Cox model confirmed that effects operated within genes rather than between them. Leave-one-gene-out analysis identified the contribution of each gene.

## Sensitivity analyses

Sensitivity analyses 1a and 2–7 were pre-specified in the original analysis plan: (1a) germline-enriched versus somatic-enriched carriers [23]; (2) restricted mean survival time at τ = 60, 81, 120 months [24]; (3) progressive tumor purity filtering (≥ 0.2 through ≥ 0.5); (4) E-value for unmeasured confounding [25]; (5) BRCA1 promoter methylation; (6) TMB and immune correlates [31]; (7) stage-adjusted Cox. The predictor benchmark, concordance analysis (1b), tier-level survival, and exploratory CCF were added during a pre-submission audit and are treated as secondary comparative analyses. An exploratory cancer cell fraction analysis using an approximate VAF/purity metric is reported in the Supplementary Appendix.

## Software and reproducibility

Analyses used Python 3.12.3 (seed = 42) with pinned package versions (requirements.txt). The pipeline is available at https://github.com/rm147747/AlphaHRD under CC BY 4.0.
