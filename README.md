# PhenoAge Multi-Layer Biological Interpretation Map — Replication Guide

Everything needed to understand and reproduce the analysis from raw public
data. Runtime environment: R 4.4.x, Python 3 (stdlib only for the demo
builder), MAGMA v1.10, ~16 GB RAM, ~10 GB disk. Total compute ≈ 6–8 h
(dominated by MAGMA gene analysis, ~25–40 min per trait single-core;
parallelizable).

---

## What this is (one paragraph)

This package contains a multi-layer biological interpretation map of PhenoAge,
the blood-test-based ageing measure (Levine 2018; Kuo 2021 constants). Layer 1
(T1, exact) decomposes the PhenoAge equation into its 10 clinical terms and
computes each term's contribution to a person's ageing gap. Layer 2 (T2,
genetic) links each biomarker to genes and pathways via MAGMA gene-set
analysis of UK Biobank biomarker GWAS (Neale Lab round 2), adds
Sinnott-Armstrong 2021 fine-mapped credible-set variants, and adds an
inverse-variance Mendelian-randomization layer testing whether
genetically-predicted biomarker levels causally influence PhenoAge
acceleration (Kuo GWAS as outcome). Layer 3 (T3, proteomic) adds 83 circulating
proteins associated with PhenoAge acceleration (Bloch 2026 preprint). Layer 4
(T4, network) adds PPI interactions and network-proximity tests between
protein sets and hallmarks of ageing. Layer 5 (T5, annotation) attaches
OpenGenes hallmark-of-ageing labels with confidence scores. Every edge is
typed, evidence-graded (T1–T5), and cited; nothing is a causal claim except
the explicitly-caveated MR layer. The primary artifact is `graph.json`
(1,147 nodes, 1,626 edges) plus flat CSVs, ready to load into a custom app;
`phenoage_map_demo.html` is a self-contained interactive visualization with a
synthetic example person (age 58, PhenoAge 66.5, acceleration +8.49 y).

---

## Package contents

| File | What it is |
|---|---|
| `graph.json` | Primary artifact: typed, cited, tier-graded evidence graph |
| `nodes.csv`, `edges.csv`, `citations.csv`* | Flat views of the graph (*citations.csv is embedded in graph.json; regenerate via assemble_graph.R) |
| `SCHEMA.md` | Every field, tier, relation, threshold, and caveat |
| `person_example_decomposition.csv`, `person_example_summary.csv`, `reference_person.csv` | T1 worked example (synthetic person) |
| `gwas_genes_<trait>.csv` | Full MAGMA gene-level results per trait (10 traits) |
| `magma_pathways_<trait>.csv` | Full MAGMA gene-set results per trait (FDR added post-hoc) |
| `mr_results.csv`, `mr_instruments.csv` | IVW MR estimates and per-instrument ratios |
| `network_proximity_bloch.csv` | Protein-set vs hallmark proximity tests |
| `bloch_diablo_features.csv`, `bloch_proteins_univariate.csv` | Parsed Bloch 2026 supplementary tables |
| `figures/` | Evidence-depth heatmap, pathway dot plot, Fig 3 replication (longevity module + 11 hallmark modules; SVG + PNG) |
| `phenoage_map_demo.html` | Self-contained interactive demo (open in any browser) |
| `scripts/` | All code, in pipeline order (below) |

---

## Data sources (all public)

1. **Neale Lab UK Biobank round 2 GWAS** — `s3://broad-ukb-sumstats/round2/additive-tsvs/{field_id}_irnt.gwas.imputed_v3.both_sexes.tsv.bgz`
   Field IDs used: WBC 30000, MCV 30040, RDW 30070, lymphocyte % 30180,
   albumin 30600, ALP 30620, creatinine 30700, CRP 30710, glucose 30740.
   (Verify IDs at biobank.ndph.ox.ac.uk — several are counterintuitive.)
   Variant IDs are `chr:pos:ref:alt`; convert to rsIDs via the 1000G EUR bim.
2. **Kuo 2021 PhenoAgeAccel GWAS sumstats** — figshare file 23739971
   (CC BY 4.0), 16.4M SNPs, N=107,460. Paper: PMC8208797 (the Wiley
   supplementary is 403-blocked; the equation is quoted in the PMC text).
3. **Sinnott-Armstrong 2021 fine-mapped associations** — public Google Sheet
   linked from github.com/rivas-lab/biomarkers (`fine_mapping/filtered_hits`);
   3,433 credible-set variants (posterior ≥ 0.99) across 35 lab biomarkers.
   (The NIH figshare mirror is WAF-blocked; the Google Sheet works.)
   **Column note: `beta` is the effect of allele2 (ALT), not allele1.**
4. **Bloch 2026** — medRxiv doi 10.64898/2026.03.02.26347435, supplementary
   PDF (Tables S5 + DIABLO feature tables). Preprint — flag as such.
5. **Longevity module** (Gross 2025) — github.com/BnayaGross/Longevity-module:
   `PPI_2022.csv` (531,597 edges) and `Gene_hallmarks.csv` (3,079
   gene-hallmark annotations with confidence 1–5).
6. **MAGMA v1.10** — ctg.cncr.nl is broken (redirects to HTML); recover
   binary + `g1000_eur.zip` + `NCBI37.3.zip` from Wayback Machine snapshots
   (20220402124539, 20190427120156, 20190427120021).
7. **MSigDB** — via the R `msigdbr` package (Hallmark + C2:CP:KEGG + C5:GO:BP;
   7,588 sets).
8. **NHANES IV reference medians** — for the reference person in the T1 layer.

---

## Pipeline (run in this order)

| Step | Script | Inputs → outputs | Notes |
|---|---|---|---|
| 0 | — | Download all sources above into `downloads/` | ~10 GB |
| 1 | `scripts/phenoage_module.R` | Kuo/Levine constants → `phenoage_module.R` functions | Implements both constant sets (kuo2021 primary, levine2018 unrounded; ~1.6 y offset). Validated: exact linearity, exact age-delta, r=1.0 vs BioAge package on NHANES IV |
| 2 | `scripts/person_example.R` | Step 1 → `person_example_*.csv`, `person_example.rds` | Synthetic person; contributions sum exactly to the ageing gap |
| 3 | `scripts/parse_bloch.R` | Bloch supp PDF (needs `pdftotext -layout`, from poppler-utils) → `bloch_*.csv` | pdftools R package fails to install; use the CLI |
| 4 | `scripts/make_genesets.R` | msigdbr → `magma/genesets_magma.txt` | 7,588 sets, column format `set<TAB>ncbi_gene` |
| 5 | `scripts/run_magma.sh` (or `run_magma_parallel.sh`) | Neale bgz + Kuo txt → per-trait `_pvalues.txt`, `.genes.out`, `.genes.raw`, `_gs.gsa.out` | See gotchas below |
| 6 | `scripts/integration_network.R` | Bloch proteins + PPI → `network_proximity_bloch.csv`, `ppi_graph.rds` | Multi-source BFS mean-min distance vs 200 degree-binned permutations |
| 7 | `scripts/run_mr.R` | Fine-mapped hits + Kuo sumstats → `mr/mr_results.csv`, `mr/mr_instruments.csv` | One instrument per region; IVW; harmonization verified vs APOE rs429358 |
| 8 | `scripts/assemble_graph.R` | All above → `graph.json`, `nodes.csv`, `edges.csv` | Thresholds: genes p<1e-6 (top 50/trait), pathways FDR<0.05 (top 30/trait) |
| 9 | `scripts/export_magma_tables.R` | MAGMA outputs → `gwas_genes_*.csv`, `magma_pathways_*.csv` | |
| 10 | `scripts/make_figures.R` | Graph + MAGMA outputs → `figures/*` | ggplot2 + svglite |
| 11 | `build_demo.py` (in session workspace) | graph.json + cytoscape.min.js → `phenoage_map_demo.html` | Embeds cytoscape.js (cdnjs) + graph into one file |
| 12 | `scripts/fig3_longevity_module.R` | PPI + hallmark annotations → `figures/fig3l_longevity_module.*`, `figures/fig3ak_hallmark_modules.*` | Replicates Gross et al. 2026 Fig 3: per-hallmark LCCs agglomerated into the longevity module; node color = hallmark (black = multi-hallmark), size = n hallmarks, labels ≥ 5 hallmarks. Note: paper's Fig 3 panels are smaller (359 genes total) than the full-set LCCs (~1,088); their exact gene subset is not specified in the main text |

### MAGMA gotchas (each cost us a run)

- `--gene-annot` takes the **annotation output** of a prior
  `--annotate window=10,1.5 --snp-loc <bim> --gene-loc <gene.loc>` step —
  not the gene.loc itself.
- The window modifier belongs on `--annotate`, not `--gene-loc`.
- p-value input: `--pval use=1,2 N=<n>` (capital N; no `col=` syntax).
- Gene-set input is row-based by default; for a two-column
  `set<TAB>gene` file add `col=2,1` to `--set-annot`.
- `awk` cannot read bgzf directly — pipe through `zcat`.
- Filter NaN p-values before building pfiles.
- Gene-set output is `.gsa.out` (not `.genes.out`); FDR is not included —
  compute BH post-hoc.

### MR gotchas

- The Sinnott-Armstrong `beta` column is the **ALT-allele (allele2) effect**.
  Verify harmonization against a known-direction locus (we used APOE
  rs429358: e4 allele must give positive PhenoAgeAccel).
- Both GWAS are UKB-based — sample overlap may inflate estimates; Q
  heterogeneity is significant for 4/5 biomarkers. Treat as supportive,
  not definitive.

---

## Validation checks (all pass)

- PhenoAge module: exact linearity, exact age-delta additivity, r=1.0 and
  median |error| 0.001 y vs the BioAge package on NHANES IV.
- Graph: 0 dangling edges, 0 edges without citation, tiers ⊆ {T1..T5}.
- MR: sign verification vs APOE rs429358; beta correlations per trait all
  biologically directed after ALT-allele correction.
- Demo HTML: embeds the exact graph.json (node/edge counts verified).

## Known limitations

1. Non-causal by construction (MR layer is the only causal-inference content
   and carries its own caveats).
2. Bloch 2026 is a non-peer-reviewed preprint.
3. CRP unit ambiguity (mg/dL vs mg/L) between Bloch S5 and the Kuo equation.
4. Two PhenoAge constant sets exist (~1.6 y offset); kuo2021 used throughout
   because the composite GWAS used it.
5. European-ancestry LD reference only.
6. Sinnott-Armstrong fine-mapping/MR cover the 5 biochemistry biomarkers
   only (their 35 traits exclude blood counts).
7. The example person is synthetic.

## References

**Core measures**

1. Liu Z, Kuo P-L, Horvath S, Crimmins E, Ferrucci L, Levine ME (2018). A new aging measure captures morbidity and mortality risk across diverse subpopulations from NHANES IV: A cohort study. *PLoS Med* 15(12):e1002718. <https://doi.org/10.1371/journal.pmed.1002718>
2. Levine ME, Lu AT, Quach A, et al. (2018). An epigenetic biomarker of aging for lifespan and healthspan. *Aging* 10(4):573–591. <https://doi.org/10.18632/aging.101414>

**Molecular evidence layers**

3. Bloch S, de Magalhães JP (2026). Proteomic and metabolomic signatures of PhenoAge acceleration (UK Biobank, DIABLO multi-omics). Preprint. <https://www.medrxiv.org/content/10.64898/2026.03.02.26347435v1>
4. Kuo P-L, Horvath S, Crimmins EM, Levine ME (2021). Genetic associations for two biological age measures point to distinct aging phenotypes. *Aging Cell* 20(6):e13376. <https://doi.org/10.1111/acel.13376> — GWAS summary statistics: [PhenoAgeAccel](https://doi.org/10.6084/m9.figshare.12620291.v1), [BioAgeAccel](https://doi.org/10.6084/m9.figshare.12620366.v1)
5. Sinnott-Armstrong N, Tanigawa Y, Amar D, et al. (2021). Genetics of 35 blood and urine biomarkers in the UK Biobank. *Nat Genet* 53:185–194. <https://doi.org/10.1038/s41588-020-00757-z>
6. Gross B, Ehlert J, Gladyshev V, Loscalzo J, Barabási A-L, et al. (2026). Network-driven discovery of repurposable drugs targeting hallmarks of aging. *Nature Aging*. <https://www.nature.com/articles/s43587-026-01161-8>

**Resources & methods**

7. López-Otín C, Blasco MA, Partridge L, Serrano M, Kroemer G (2023). Hallmarks of aging: An expanding universe. *Cell* 186(1):243–278. <https://doi.org/10.1016/j.cell.2022.12.045>
8. Rafikova E, et al. (2024). Open Genes database. *Nucleic Acids Res* 52(D1):D950–D962. <https://doi.org/10.1093/nar/gkad712> — <https://open-genes.org/>
9. de Leeuw CA, Mooij JM, Heskes T, Posthuma D (2015). MAGMA: generalized gene-set analysis of GWAS data. *PLoS Comput Biol* 11(4):e1004504. <https://doi.org/10.1371/journal.pcbi.1004504>
10. Subramanian A, et al. (2005). Gene set enrichment analysis: a knowledge-based approach (MSigDB). *PNAS* 102(43):15545–15550. <https://doi.org/10.1073/pnas.0506580102>
11. Bycroft C, et al. (2018). The UK Biobank resource with deep phenotyping and genomic data. *Nature* 562:553–559. <https://doi.org/10.1038/s41586-018-0579-z>
12. Neale Lab (2018). UK Biobank GWAS. <http://www.nealelab.is/uk-biobank/>
