<h1 align="center">APHunter: Amplicon Primer Hunter</h1>
<p align="center">
  <strong>🔍 Determine the primer used in any amplicon sequencing dataset</strong>
</p>
<p align="center">
  <a href="https://anaconda.org/westraingroup/aphunter">
    <img src="https://img.shields.io/conda/v/westraingroup/aphunter.svg?style=flat-square&logo=anaconda" alt="Conda Version">
  </a>
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square" alt="License">
  </a>
  <!-- You can add more badges here, e.g., for build status, downloads, etc. -->
</p>

## Table of Contents
*   [Introduction](#introduction)
    *   [Workflow](#workflow)
        *   [Non-pool version](#non-pool-version)
        *   [Pool version](#pool-version)
*   [User Manual](#user-manual)
    *   [Installation and Execution](#installation-and-execution)
    *   [Parameters](#parameters)
    *   [Pool or Non-pool?](#pool-or-non-pool)
    *   [Output Files](#output-files)
    *   [Understanding ITS Primer 'Siblings'](#understanding-its-primer-siblings)
    *   [Verifying Results (Optional)](#verifying-results-optional)
*   [References](#references)
*   [License](#license)
*   [Citation](#citation)
*   [Contact](#contact)
*   [Future Updates](#future-updates)

## Introduction
APHunter is a Python-based tool developed for the identification of primers in amplicon sequencing data from `.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz` files. It currently offers two operational versions: a "non-pool" mode and a "pool" mode, each tailored for different analytical scenarios. A detailed comparison between these modes can be found in the [Pool or Non-pool?](#pool-or-non-pool) section.

### Workflow
The core workflow of APHunter involves several key steps, with variations depending on the chosen mode:

#### Non-pool version
1.  **Reads Truncation:** All input reads are truncated to the first 100 base pairs using SeqKit.
2.  **Quality Filtering:** Low-quality reads are removed using VSEARCH, employing a maximum expected error (maxEE) threshold of 1.0.
3.  **Sample Filtering:** Samples with fewer than 1000 reads post-filtering are discarded.
4.  **Base Frequency Analysis:** For each remaining sample, the frequency of each nucleotide (A, T, C, G) at every base position across the truncated reads is calculated and visualized.
5.  **Consensus Sequence Generation:** A consensus sequence is derived for each sample based on the most frequent base (MFB) at each position. The consensus sequence is truncated when five consecutive MFBs exhibit a frequency below 0.75, indicating a potential transition into a less conserved region.
6.  **Primer Identification:** The generated consensus sequences are aligned against a primer database (e.g., an ITS primer database) using BLASTn-short. Significant matches identify potential primers used in the amplification.

#### Pool version
1.  **Sample Pooling:** All input samples are pooled into a single file.
2.  **Reads Truncation:** All reads in the pooled file are truncated to the first 100 base pairs using SeqKit.
3.  **Quality Filtering:** Low-quality reads are removed using VSEARCH, employing a maximum expected error (maxEE) threshold of 1.0.
4.  **Reads Clustering:** The filtered reads are clustered at a likelihood threshold of 0.9 using VSEARCH.
5.  **Cluster Filtering:** Clusters containing fewer than 0.05% of the total reads, or fewer than 100 reads (if 0.05% of total reads is less than 100), are discarded.
6.  **Base Frequency Analysis:** For each remaining cluster, the frequency of each nucleotide (A, T, C, G) at every base position across the reads within that cluster is calculated and visualized.
7.  **Consensus Sequence Generation:** A consensus sequence is derived for each cluster based on the MFB at each position. The consensus sequence is truncated when five consecutive MFBs exhibit a frequency below 0.75.
8.  **Primer Identification:** The consensus sequences from clusters are aligned against a primer database (e.g., an ITS primer database) using BLASTn-short to identify potential primers.

## User Manual
This section provides guidance on installing and using APHunter.

### Installation and Execution
APHunter is available on Anaconda and can be installed via the `westraingroup` channel:
```bash
conda install westraingroup::aphunter
```
It is recommended to install APHunter in a new conda environment, as this will automatically set up all required dependencies.

Once installed, APHunter can be run with the following basic command structure:
```bash
aphunter -i /your/input/directory -o /your/output/directory -s your_sample_file_suffix
```

### Parameters
The following table details the command-line parameters available for APHunter:

| Parameter               | Description                                                                                                                                                                |
| :---------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-i`, `--input`           | **Input folder.** Directory containing `.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz` files. (Required)                                                                      |
| `-o`, `--output`          | **Output folder.** Directory to save the results. If the folder does not exist, it will be created automatically. (Required)                                             |
| `-s`, `--suffix`          | **File suffix.** Specify the common suffix of the target filenames to be processed (e.g., `_1.fastq.gz`, `_2.fastq`). (Required)                                           |
| `-p`, `--pool`            | **Pooled mode.** Enable this option to use the pooled version of APHunter. (Optional, defaults to non-pool mode)                                                        |
| `-t`, `--threads`         | **Number of threads.** Specifies the number of threads for parallel processing. Defaults to 8 or half of the system’s maximum thread count, whichever is greater. (Optional) |
| `-e`, `--e_value`         | **E-value threshold.** The E-value cutoff for filtering BLAST results. Default is `1e-3`. Lower values indicate higher stringency. (Optional)                               |
| `-h`, `--help`            | **Help.** Display usage information and descriptions of all available options.                                                                                             |

### Pool or Non-pool?
Choosing between the pool and non-pool versions depends on your specific research question and dataset characteristics.

| Feature                             | Non-pool version                                                                                   | Pool version                                                                                                                                  |
| :---------------------------------- | :------------------------------------------------------------------------------------------------- | :-------------------------------------------------------------------------------------------------------------------------------------------- |
| Workflow Link                       | [View Non-pool Workflow](#non-pool-version)                                                        | [View Pool Workflow](#pool-version)                                                                                                           |
| Output Focus                        | Information and results are sample-specific.                                                       | Information and results are cluster-specific, derived from pooled data.                                                                       |
| Execution Speed                     | Generally faster for datasets with distinct sample compositions.                                   | Can be slower due to the comprehensive clustering step (please be patient 🥺).                                                              |
| Interpretation of BLAST results     | Primer matches can be directly traced back to the specific samples in which they were detected.    | Primer matches are traced to read clusters and their relative abundance, indicating primers present across the pooled dataset.                |
| Primary Task                        | Ideal for identifying the most probable primer used in each individual sample.                     | Suited for identifying the range of primers present in a batch of samples, including detection of multiple primers if used (intentionally or not). |
| **Recommended Use Case**            | Precise, sample-level primer tracing; useful for detailed analysis and contamination tracking.     | High-throughput or batch sample analysis; efficient identification of primers across diverse clusters.                                      |

*This comparison is a guideline. We encourage users to experiment with both versions on their data and share their findings [with us](#contact) to help us further improve APHunter.*

### Output Files
APHunter generates several output files organized into subdirectories within the specified output folder.

| Directory           | File Name Pattern                    | Description                                                                                                       |
| :------------------ | :----------------------------------- | :---------------------------------------------------------------------------------------------------------------- |
| `/logs`             | `pipeline.log`                       | (Common to both versions) Comprehensive log file detailing all processing steps and messages.                     |
| `/filtered_samples` | `{sample_name}_filtered.fasta`       | (Non-pool version only) FASTA file of filtered reads for each sample, ready for subsequent analysis.              |
| `/stats`            | `{sample_name}_base_freqs.csv`       | (Non-pool version only) CSV file containing base frequency information for each filtered sample.                  |
| `/plots`            | `{sample_name}_freq_plot.png`        | (Non-pool version only) Line plot visualizing base frequencies for each filtered sample.                          |
| `/blast`            | `sample_consensus_for_blast.fasta`   | (Non-pool version only) FASTA file of consensus sequences derived from each sample, used as BLAST input.          |
| `/blast`            | `per_sample_blast_output.tsv`        | (Non-pool version only) Raw BLASTn output in tabular format for per-sample consensus sequences.                   |
| `/blast`            | `per_sample_final_blast_results.csv` | (Non-pool version only) Processed and more readable BLAST results for samples.                                    |
| `main output dir`   | `all_centroids.fasta`                | (Pool version only) FASTA file containing all representative centroid sequences identified by VSEARCH clustering. |
| `/clusters`         | `{cluster_name}.fasta`               | (Pool version only) Individual FASTA files containing reads belonging to each filtered cluster.                   |
| `/stats`            | `{cluster_name}_base_freqs.csv`      | (Pool version only) CSV file containing base frequency information for each cluster.                              |
| `/plots`            | `{cluster_name}_freq_plot.png`       | (Pool version only) Line plot visualizing base frequencies for each cluster.                                      |
| `/blast`            | `cluster_consensus_for_blast.fasta`  | (Pool version only) FASTA file of consensus sequences derived from each cluster, used as BLAST input.             |
| `/blast`            | `pooled_blast_output.tsv`            | (Pool version only) Raw BLASTn output in tabular format for pooled cluster consensus sequences.                   |
| `/blast`            | `pooled_final_blast_results.csv`     | (Pool version only) Processed and more readable BLAST results for clusters.                                       |

*Note: `{sample_name}` refers to the original input file prefix, and `{cluster_name}` (e.g., `cluster.1`, `cluster.2`) refers to the identifier assigned to each cluster.*

### Understanding ITS Primer 'Siblings'
Due to sequence similarities, certain primers within the ITS primer database (or any comprehensive primer database) may be co-detected by BLAST. We term these "primer siblings." APHunter provides a pre-compiled list of such ITS primer groups to aid in result interpretation.

These groups were identified by an all-versus-all BLAST search within the ITS primer database. An alignment pair was considered significant if the alignment length was ≥15 base pairs with a percent identity >90% (resulting in an E-value < 1e-3).

**NOTE:** Primers are grouped if any primer in the group shares sufficient similarity with at least one other member of the same group. This means not all pairs within a group are necessarily highly similar to each other directly, but they are connected through a chain of similarity.

*   ITS1_1_fwd - ITS1ngs_1_fwd - ITS1Fngs_1_fwd - ITS1F_1_fwd - ITSOF-T_1_fwd - ITS1F_KYO1_1_fwd
*   58A1F_2_fwd - 58A2F_2_fwd - ITS3_2_fwd - 58A2R_1_rev - ITS2_1_rev
*   fITS7_2_fwd - gITS7_2_fwd - ITS86F_2_fwd - fITS7r_1_rev
*   ITS4_2_rev - ITS4-Fun_2_rev - ITS4ngs_2_rev
*   ITS9mun_1_fwd - ITS9MUNngs_1_fwd
*   ITS4_KYO2_2_rev - ITS4-Tul2_2_rev
*   ITS1-F_KYO2_1_fwd - ITS5_1_fwd
*   ITS4-Clav_2_rev - LR0B_2_rev
*   BITS_1_fwd - ITS1catta_1_fwd
*   ITS4B1_2_rev - ITS4B_2_rev
*   LR3_2_rev - TW13_2_rev

**Reminder:** When interpreting the results, there are four main senarios. You may refer to the following table to quickly identify the corresponding type of primer possibly used in your dataset:

| ITS region\Read orientation | Fwd (`_1` or `.1` input) | Rev (`_2` or `.2` input) |
| --------------------------- | ------------------------ | ------------------------ |
| **ITS1**                    | `_1_fwd` primers         | `_1_rev` primers         |
| **ITS2**                    | `_2_fwd` primers         | `_2_rev` primers         |

### Verifying Results (Optional)
To further validate primers identified by APHunter, users can employ tools like [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) to quantify the proportion of reads containing a specific primer sequence. This can provide an estimate of the primer's coverage in the raw or filtered dataset.

An example Cutadapt command (removing 5' end adapter) is:
```bash
cutadapt -g PRIMER_SEQ --minimum-length 1 --discard-untrimmed -o matched.fastq input.fastq > cutadapt.log
```
Put your 'suspect' primer into `PRIMER_SEQ`, adjust the error tolerance (`-e`, default is 0.1), input (`input.fastq`), and output (`-o matched.fastq`) parameters as needed.

## References
APHunter leverages several established bioinformatics tools and curated databases.

**External Tools:**
*   [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
*   [SeqKit](https://bioinf.shenwei.me/seqkit/)
*   [VSEARCH](https://github.com/torognes/vsearch)

**Database Sources (primarily for ITS primers):**
1.  Bokulich, N. A., & Mills, D. A. (2013). Improved selection of internal transcribed spacer-specific primers enables quantitative, ultra-high-throughput profiling of fungal communities. *Applied and Environmental Microbiology*, *79*(8), 2519–2526. https://doi.org/10.1128/AEM.03870-12
2.  Nilsson, R. H., Anslan, S., Bahram, M., et al. (2019). Mycobiome diversity: high-throughput sequencing and identification of fungi. *Nature Reviews Microbiology*, *17*, 95–109. https://doi.org/10.1038/s41579-018-0116-y
3.  UNITE Community. (n.d.). Primers for amplifying and sequencing fungal rDNA ITS region. Retrieved from https://unite.ut.ee/primers.php

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Citation
If you use APHunter in your research, please cite it as follows (APA style):

Zhou, W., Wang, X., & Zheng, J. S. (2025). *APHunter* (Version 1.0) [Computer software]. Westlake University. https://github.com/WeStrainGroup/APHunter

## Contact
For questions, suggestions, or bug reports, please contact us via email:
*   Wenhao Zhou: [zhouwenhao@westlake.edu.cn](mailto:zhouwenhao@westlake.edu.cn)
*   Xinyu Wang: [wangxinyu30@westlake.edu.cn](mailto:wangxinyu30@westlake.edu.cn)

We greatly appreciate your feedback and contributions!

## Future Updates
We plan the following enhancements for APHunter (based on current version 1.1.0):
1.  **Improved Log Output (Target: v1.1.1):** More structured and informative logging.
2.  **BLAST Result Prioritization (Target: v1.2.0):** Enhanced analysis of BLAST results to provide a ranked recommendation of the "top suspect" primer(s).
3.  **User-Supplied BLAST Database (Target: v1.3.0):** Functionality to allow users to provide their own custom primer databases for BLAST.

If you have ideas regarding these or other potential features, please do not hesitate to [contact us](#contact).
