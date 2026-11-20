# PURPLE-seeded AmpliconSuite

Nextflow workflow for GRCh38 paired-end tumour WGS BAMs with matching indexes and
PURPLE results. The container pins AmpliconSuite 1.6.2, AmpliconArchitect 1.6.r2
and AmpliconClassifier 2.0.1. AA uses Clarabel.

## Modes

`--mode purple` (default) seeds AA from PURPLE as described below.

`--mode default` runs AmpliconSuite with its own default seeding: CNVkit on the
tumour BAM with the matched normal BAM (`--normal_bam`), then the wrapper's
`cnv_prefilter` and `amplified_intervals.py` filters. No other AmpliconSuite or
AA option is changed from its default. The input CSV needs
`sample,bam,bai,normal_bam,normal_bai`; `purple` and `mask` are not used.

| Parameter | Input |
| --- | --- |
| `mode` | `purple` (default) or `default`. |
| `input` | CSV with unique `sample,bam,bai,kras_qc,kras_vaf,kras_depth,kras_variant` columns (purple mode), or `sample,bam,bai,normal_bam,normal_bai` (default mode). |
| `purple` | ZIP with `<sample>.purple.purity.tsv` and `<sample>.purple.cnv.somatic.tsv` at its root. |
| `mask` | `GRCh38_merged_centromeres_conserved_sorted.bed` from the same AA reference archive. |
| `reference` | Official AA GRCh38 tarball with a top-level `GRCh38/` directory. |
| `outdir` | Writable output directory or S3 prefix. |

Use a manifest reconciled against the BAM and PURPLE provenance. The KRAS fields
carry optional precomputed QC and can be empty. Check fit and sequencing QC
before including samples in the manifest.

PURPLE absolute tumour copy number divided by sample ploidy defines gain
(`1.5 ≤ ratio < 3`) and amplification (`ratio ≥ 3`). Autosomal and sex-chromosome
coordinates are converted from PURPLE's 1-based inclusive format to BED. The
workflow reports all gains and amplifications. Amplifications are seeded at
any length. Gains are seeded only when tumour CN is at least 4.5 and the
contiguous block of such gains is shorter than 10 Mb; gain blocks are assessed
separately from amplified blocks. The AA conserved-gain and centromere mask is
subtracted from every seed. It passes the remaining intervals in a
`*_AA_CNV_SEEDS.bed` file that bypasses CNVKit and the additional default
AmpliconSuite seed filter. AA and AC still use their native reconstruction and
classification settings. The fitted tumour CN sets seed eligibility. AA
estimates graph copy numbers from the underlying reads.

The output under `seed_audit/` includes per-sample calls, exact seed intervals
and sample status, including `no_eligible_seed`. Seeded samples have an AA
summary, AA/AC outputs, and a status file under `samples/<sample>/`. The
feature-level origin table retains the native ecContext label and groups
the labels according to the AmpliconSuite paper. Origin groups are
structure-based predictions. Use the original ecContext labels when reviewing
two-foldback and BFB-like events. No eligible seed does not establish ecDNA
absence.
