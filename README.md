<h1 align="center">APHunter: Amplicon Primer Hunter</h1>
<p align="center">
  <strong>Determine the primer used in any amplicon sequencing dataset.</strong>
</p>


## Introduction

APHunter identifies which reference primer best explains the conserved 5′ end of amplicon reads. It accepts `.fastq`, `.fq`, `.fastq.gz`, and `.fq.gz` files and supports two analysis modes:

- **Per-sample (default):** build and query one consensus per sample.
- **Pooled (`-p`):** merge all samples, cluster the reads, and query one consensus per retained cluster.

## Installation

APHunter is distributed and supported through Conda. Install it with:

```bash
conda create -n aphunter westraingroup::aphunter
conda activate aphunter
```

The Conda package installs the Python and command-line dependencies and configures the bundled ITS database as the default. Both built-in databases remain selectable with `--blast-db`.

## Command-line options

| Option | Description |
| ------ | ----------- |
| `-i`, `--input` | Input directory containing FASTQ files. Required. |
| `-o`, `--output` | Output directory. Created when necessary. Required. |
| `-s`, `--suffix` | Filename suffix used to select inputs and infer read orientation. Required. |
| `-p`, `--pool` | Enable pooled clustering mode. |
| `-t`, `--threads` | Worker and BLAST threads. Default: `max(8, floor(CPU count / 2))`. |
| `-e`, `--e_value` | BLAST E-value ceiling. Default: `1e-3`. |
| `--blast-db {its,16s}` | Select a bundled primer database. |
| `--blast-db-dir DIR` | Use a custom nucleotide BLAST index directory. |
| `-h`, `--help` | Show command-line help. |

Example using the bundled ITS database in pooled mode:

```bash
aphunter -i /path/to/fastq_dir -o /path/to/output -s _R1.fastq.gz --blast-db its -p
```

## How it works

1. Select input files whose names end with the requested suffix.
2. Truncate reads to positions 1–100 with SeqKit.
3. Filter reads with VSEARCH at `maxEE = 1.0`.
4. In pooled mode, cluster reads at 90% identity and retain clusters with at least `max(100, N/2000)` reads, where `N` is the pooled filtered-read count. In per-sample mode, discard samples with fewer than 1,000 filtered reads.
5. Count per-position base frequencies with weighted IUPAC ambiguity support.
6. Build a consensus using a 75% dominant-base threshold. Isolated low-confidence positions may be retained, but consensus extension stops when low-confidence signal persists across the following two positions.
7. Plot base frequencies and compare non-empty consensus sequences against the selected primer database with BLAST+.

## Output layout

Numbered directories contain intermediate results; `final/` contains the primary deliverables.

| Mode | Intermediate flow |
| ---- | ----------------- |
| Per-sample | `1_filtered_samples/` → `2_stats/` → `3_plots/` → `4_blast_work/` |
| Pooled | `1_pool_work/` → `2_clusters/` → `3_stats/` → `4_plots/` → `5_blast_work/` |

Both modes also create:

- `final/*_final_blast_results.csv`: ranked query–hit rows.
- `final/*_top_suspect.csv`: one best-hit or miss row per query.
- `final/primer_call_summary.txt`: run-level likely primer, support, identity, confidence, and per-primer breakdown.
- `logs/pipeline.log`: complete execution log.

Base-frequency CSV files report `Pos`, `MaxBase`, `MaxFreq`, `InformativeReads`, `AmbiguousReads`, `NReads`, and `UnrecognizedReads`.

## BLAST databases and matching

| Selection | Content |
| --------- | ------- |
| Default or `--blast-db its` | Fungal ITS primers |
| `--blast-db 16s` | Bacterial 16S primers |

The bundled ITS database is used when no database option is given. Use `--blast-db 16s` for the bundled bacterial library or `--blast-db-dir DIR` for a custom directory containing one nucleotide BLAST index. Conda configures the bundled database location automatically.

BLAST uses `strand=both`, `word_size=7`, and the requested E-value ceiling. Queries shorter than 50 bp use `blastn-short`; longer queries use `blastn`.

After the initial BLAST ranking, APHunter uses recognized R1/R2 markers in `-s` to infer read orientation. It removes hits explicitly labeled as the opposite orientation and then reranks the remaining hits. For 16S IDs, direction is identified by a terminal `f` or `r`; for ITS IDs, it is identified by `_fwd` or `_rev`. Hits whose primer direction cannot be classified are retained. If read orientation cannot be inferred from `-s`, orientation filtering is skipped.

## Primer identifiers and APVisual

APHunter reports primer database identifiers in its result tables and summary. Look up an identifier in [APVisual](https://apvisual.netlify.app/) to view its sequence, position, citation, and reference link.

### 16S identifiers

For bacterial 16S primers, APHunter uses coordinate-normalized identifiers such as `335f` and `786r`. These identifiers correspond to the labels shown in APVisual and encode the curated reference position and orientation. Renaming changes only the identifier, not the primer sequence.

### ITS identifiers

ITS subject IDs include the target region and direction. In suffixes such as `_1_fwd` and `_2_rev`, `1` and `2` indicate ITS1 and ITS2, while `fwd` and `rev` indicate primer orientation. APVisual displays the corresponding base name without this suffix.

## Choosing an analysis mode

| Consideration | Per-sample | Pooled |
| ------------- | ---------- | ------ |
| Analysis unit | Individual sample | Read cluster across all samples |
| Best suited for | Sample-level tracing and contamination checks | Batch-level detection of one or more primer signals |

Runtime depends on the number of files, total reads, and sequence diversity; pooled mode adds a global clustering step.

## Interpreting similar ITS primers

Closely related ITS primers can produce similar BLAST hits. The groups below were derived from an all-versus-all search in which links required an alignment length of at least 15 bp and identity above 90%. Groups are connected components, so not every pair within a group must match directly. APHunter reports individual subjects and does not merge these groups automatically.

<details>
<summary>Show ITS primer sibling groups</summary>

- 58A1F_2_fwd – 58A2F_2_fwd – ITS3_2_fwd
- 58A2R_1_rev – ITS2_1_rev
- BITS_1_fwd – ITS1catta_1_fwd
- fITS7_2_fwd – gITS7_2_fwd – ITS86F_2_fwd
- ITS1_1_fwd – ITS1ngs_1_fwd
- ITSF_1_fwd – ITS1F_KYO1_1_fwd – ITSFngs_1_fwd – ITSOF_1_fwd
- ITS1F_KYO2_1_fwd – ITS5_1_fwd
- ITS4_2_rev – ITS4ngs_2_rev – ITS4-Fun_2_rev
- ITS4_KYO2_2_rev – ITS4-Tul2_2_rev
- ITS4-Clav_2_rev – LR0B_2_rev
- ITS4B_2_rev – ITS4B1_2_rev
- ITS4PL_2_rev – LR21-Ath_2_rev
- ITS9mun_1_fwd – ITS9MUNngs_1_fwd
- LR21_2_rev – NL6Bmun_2_rev
- LR3_2_rev – TW13_2_rev

</details>


## Optional validation

[Cutadapt](https://cutadapt.readthedocs.io/en/stable/) can estimate how many reads contain a suspected primer at the 5′ end. Cutadapt is optional and is not installed with APHunter:

```bash
cutadapt -g PRIMER_SEQ --minimum-length 1 --discard-untrimmed -o matched.fastq input.fastq > cutadapt.log
```

Replace `PRIMER_SEQ` and adjust the error tolerance (`-e`, default `0.1`) as needed.

## References

External tools: [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [SeqKit](https://bioinf.shenwei.me/seqkit/), and [VSEARCH](https://github.com/torognes/vsearch).

Primer resources: see [APVisual](https://apvisual.netlify.app/) for the primer-to-publication references associated with the curated ITS and 16S libraries.

## Citation

Zhou, W., Wang, X., & Zheng, J. S. (2025). *APHunter* (Version 2.0.0) [Computer software]. Westlake University. [https://github.com/WeStrainGroup/APHunter](https://github.com/WeStrainGroup/APHunter)

## License and contact

APHunter is distributed under the [CC BY-NC-SA 4.0 license](LICENSE).

For questions, bug reports, or primer updates, contact Wenhao Zhou ([zhouwenhao@westlake.edu.cn](mailto:zhouwenhao@westlake.edu.cn)) or Xinyu Wang ([wangxinyu30@westlake.edu.cn](mailto:wangxinyu30@westlake.edu.cn)).
