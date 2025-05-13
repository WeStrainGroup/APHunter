#!/usr/bin/env python
# -- coding: utf-8 --

import os
import subprocess
import collections
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import multiprocessing
import time
import datetime
import gzip
import argparse
import glob
import functools
import sys

# ====================== Parameter Settings ======================
THREADS = max(8, multiprocessing.cpu_count() // 2 if multiprocessing.cpu_count() else 16)
SEQKIT_RANGE = "1:100"
MIN_READS_FLOOR = 100
MAX_EE = 1.0
MIN_READS_PER_SAMPLE_NON_POOL = 1000
LOG_INCREMENT_LOOPS = 25
LOG_INCREMENT_PARALLEL = 50
CLUSTER_IDENTITY_THRESHOLD = 0.9
BLAST_EVALUE = 1e-3
BLAST_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore"
BLAST_HEADER_MAP = {
    "qseqid": "query_id",
    "sseqid": "subject_id",
    "pident": "%_identity",
    "length": "alignment_length",
    "mismatch": "mismatches",
    "gapopen": "gap_opens",
    "qlen": "query_length",
    "qstart": "q_start",
    "qend": "q_end",
    "slen": "subject_length",
    "sstart": "s_start",
    "send": "s_end",
    "evalue": "evalue",
    "bitscore": "bit_score"
}
PLOT_DPI = 300
COLOR_MAP = {'A': '#FF6347', 'T': '#4682B4', 'C': '#008000', 'G': '#800080'}
THRESHOLD_FREQ = 0.75
CONSENSUS_WARNING_CONSOLE = False

# ====================== Helper Functions ======================

def setup_logging(log_file_path):
    log_dir = os.path.dirname(log_file_path)
    try:
        if log_dir: os.makedirs(log_dir, exist_ok=True)
        with open(log_file_path, "w") as f:
            f.write(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Log file initialized.\n")
        return True
    except OSError as e:
        print(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] CRITICAL Error creating log directory or initializing log file {log_file_path}: {e}. Exiting.")
        return False

def log(message, log_file_path, console_output=True):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_message = f"[{timestamp}] {message}"
    if console_output:
        print(log_message, file=sys.stdout)
        sys.stdout.flush()
    try:
        log_dir = os.path.dirname(log_file_path)
        if log_dir: os.makedirs(log_dir, exist_ok=True)
        with open(log_file_path, "a") as f: f.write(log_message + "\n")
    except Exception as e:
        err_msg = f"[{timestamp}] CRITICAL: Failed to write to log file {log_file_path}: {e}"
        print(err_msg, file=sys.stderr)
        sys.stderr.flush()

def run_command(command, log_file_path, check=True):
    log(f"Running command: {command}", log_file_path, console_output=False)
    try:
        result = subprocess.run(['/bin/bash', '-c', command], check=check, capture_output=True, text=True, errors='ignore')
        if result.stdout: log(f"Stdout: {result.stdout.strip()}", log_file_path, console_output=False)
        if result.stderr: log(f"Stderr: {result.stderr.strip()}", log_file_path, console_output=False)
        return result
    except subprocess.CalledProcessError as e:
        log(f"ERROR running command (exit code {e.returncode}): {command}", log_file_path)
        if e.stderr: log(f"Stderr: {e.stderr.strip()}", log_file_path)
        return None
    except FileNotFoundError as e:
        log(f"ERROR: Command or executable not found: {e}. Command: {command}", log_file_path)
        return None
    except Exception as e:
        log(f"ERROR: An unexpected error occurred while trying to run command: {command}\nException: {e}", log_file_path)
        return None

def count_sequences_fasta(fasta_path, log_file_path):
    if not os.path.exists(fasta_path) or os.path.getsize(fasta_path) == 0:
        return 0
    count_cmd = f"seqkit fx2tab --name --quiet '{fasta_path}' | wc -l"
    try:
        result = subprocess.run(['/bin/bash', '-c', count_cmd], check=True, capture_output=True, text=True)
        return int(result.stdout.strip())
    except (subprocess.CalledProcessError, ValueError, FileNotFoundError):
        grep_cmd = f"zgrep -c '^>' '{fasta_path}'" if fasta_path.endswith(".gz") else f"grep -c '^>' '{fasta_path}'"
        try:
            result = subprocess.run(['/bin/bash', '-c', grep_cmd], check=True, capture_output=True, text=True)
            return int(result.stdout.strip())
        except (subprocess.CalledProcessError, ValueError, FileNotFoundError):
            log(f"Warning: Fallback grep count also failed for {fasta_path}. Returning 0.", log_file_path, console_output=False)
            return 0

def count_base_frequencies(fasta_file, log_file_path):
    seqs = []
    max_len = 0
    sequence_count = 0
    valid_bases = frozenset(['A', 'T', 'C', 'G'])
    open_func = gzip.open if fasta_file.endswith(".gz") else open
    read_mode = 'rt'
    analysis_length = 0

    try:
        with open_func(fasta_file, read_mode, errors='ignore') as f:
            current_seq_parts = []
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    sequence_count += 1
                    if current_seq_parts:
                        seq_str = "".join(current_seq_parts).upper()
                        seq_str_filtered = ''.join(c for c in seq_str if c in valid_bases)
                        if seq_str_filtered:
                           seqs.append(seq_str_filtered)
                           max_len = max(max_len, len(seq_str_filtered))
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line)
            if current_seq_parts:
                seq_str = "".join(current_seq_parts).upper()
                seq_str_filtered = ''.join(c for c in seq_str if c in valid_bases)
                if seq_str_filtered:
                    seqs.append(seq_str_filtered)
                    max_len = max(max_len, len(seq_str_filtered))
    except FileNotFoundError:
        log(f"Error: Fasta file not found at {fasta_file} for base frequency counting.", log_file_path)
        return None, 0
    except Exception as e:
        log(f"Error reading base frequencies from FASTA {fasta_file}: {e}", log_file_path)
        return None, 0

    if not seqs or max_len == 0:
         return None, sequence_count
    
    analysis_length = max_len
    if ':' in SEQKIT_RANGE:
        try:
            limit_str = SEQKIT_RANGE.split(':')[-1]
            if limit_str: analysis_length = min(max_len, int(limit_str))
        except: pass

    if analysis_length <= 0:
         return None, sequence_count

    position_counters = [collections.Counter() for _ in range(analysis_length)]
    for seq in seqs:
        effective_len = min(len(seq), analysis_length)
        for i in range(effective_len):
            base = seq[i]
            position_counters[i][base] += 1
    
    freq_data = []
    for i in range(analysis_length):
        counter = position_counters[i]
        total_valid_bases_at_pos = sum(counter.values())
        if total_valid_bases_at_pos > 0:
            max_base, max_count = counter.most_common(1)[0]
            max_freq_percent = (max_count / total_valid_bases_at_pos) * 100.0
        else: max_base = 'N'; max_freq_percent = 0.0
        freq_data.append({'Pos': i + 1, 'MaxBase': max_base, 'MaxFreq': max_freq_percent})

    if not freq_data: return None, sequence_count
    return pd.DataFrame(freq_data), sequence_count

def get_cluster_id_from_path(filepath, log_file_path_for_warning=None):
    basename = os.path.basename(filepath)
    name_without_ext = basename.replace(".fasta", "").replace(".fa", "").replace(".gz", "")
    parts = name_without_ext.split('.')
    if len(parts) > 1 and parts[-1].isdigit():
        return f"Cluster_{parts[-1]}"
    else:
        clean_name = name_without_ext.replace('.', '_')
        return f"Cluster_{clean_name}"

def get_sample_id_from_path(filepath, sample_suffix_to_remove):
    basename = os.path.basename(filepath)
    sample_id = basename.replace(sample_suffix_to_remove, "")
    for ext in [".fastq.gz", ".fastq", ".fq.gz", ".fq", ".fasta.gz", ".fasta", ".fa.gz", ".fa"]:
        if sample_id.endswith(ext):
            sample_id = sample_id[:-len(ext)]
            break 
    sample_id = sample_id.replace(".", "_").replace("-", "_")
    return sample_id if sample_id else f"unknown_sample_{os.path.splitext(basename)[0]}"

def find_consensus_sequence(df, entity_id, log_file_path, freq_threshold=THRESHOLD_FREQ):
    if df is None or df.empty:
        return ""

    consensus_seq_parts = []
    threshold_percent = freq_threshold * 100.0
    lookahead_window = 2 
    started_consensus = False

    for i, row in df.iterrows():
        base, freq, pos = row["MaxBase"], row["MaxFreq"], row["Pos"]

        if base == 'N':
            if started_consensus:
                terminate = False
                if i + 1 < len(df):
                    next_freqs = df["MaxFreq"].iloc[i + 1 : min(i + 1 + lookahead_window, len(df))]
                    if next_freqs.empty or (next_freqs < threshold_percent).all():
                        terminate = True
                else:
                    terminate = True
                if terminate: break
            continue

        if not started_consensus:
            if freq >= threshold_percent:
                started_consensus = True
                consensus_seq_parts.append(base)
        else:
            if freq >= threshold_percent:
                consensus_seq_parts.append(base)
            else:
                terminate = False
                if i + 1 < len(df):
                    next_freqs = df["MaxFreq"].iloc[i + 1 : min(i + 1 + lookahead_window, len(df))]
                    if next_freqs.empty or (next_freqs < threshold_percent).all():
                        terminate = True
                    else:
                        consensus_seq_parts.append(base)
                else:
                    terminate = True
                if terminate: break
    
    final_consensus = "".join(consensus_seq_parts)
    if not final_consensus and not df.empty and CONSENSUS_WARNING_CONSOLE:
        log(f"Info: No consensus sequence generated for '{entity_id}' (no start/all below threshold).", log_file_path, console_output=True)
    return final_consensus

def plot_base_frequencies(df, entity_id, reads_in_entity, plots_dir, log_file_path):
    if df is None or df.empty:
        return False
    fig = None
    try:
        fig, ax = plt.subplots(figsize=(5, 3), dpi=PLOT_DPI)
        plot_filename = os.path.join(plots_dir, f"{entity_id}_freq_plot.png")

        ax.plot(df['Pos'], df['MaxFreq'], color='black', linewidth=0.5)
        for base, color in COLOR_MAP.items():
            subset = df[df['MaxBase'] == base]
            if not subset.empty:
                ax.scatter(subset['Pos'], subset['MaxFreq'], color=color, s=8, label=base, alpha=0.9)
        
        ax.axhline(THRESHOLD_FREQ * 100, color='red', linestyle='--', linewidth=0.6)
        ax.set_ylim(0, 125)
        plot_xmax = max(110, df['Pos'].max() + 10 if not df.empty else 0)
        ax.set_xlim(0, plot_xmax)
        ax.set_yticks(range(0, 126, 25))
        ax.set_xticks(range(0, int(plot_xmax) +1, 20 if plot_xmax > 40 else 10))
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.set_xlabel("Position", fontsize=8)
        ax.set_ylabel("Max Base Frequency (%)", fontsize=8)
        ax.set_title(f"{entity_id} | {reads_in_entity:,} reads", fontsize=9, pad=10, fontweight='bold')
        
        handles, labels = ax.get_legend_handles_labels()
        if handles:
             ax.legend(handles, labels, loc='lower right', fontsize=6, frameon=False, 
                       handletextpad=0.2, labelspacing=0.3, borderpad=0.2)
        
        plt.tight_layout(pad=0.5)
        os.makedirs(plots_dir, exist_ok=True)
        plt.savefig(plot_filename, dpi=PLOT_DPI)
        return True
    except Exception as e:
        log(f"Plotting for '{entity_id}': ERROR generating plot: {e}", log_file_path)
        return False
    finally:
        if fig: plt.close(fig)

# ====================== Core Entity Processing Functions ======================

def process_cluster_file(cluster_fasta_path, stats_dir, log_file_path):
    cluster_id = get_cluster_id_from_path(cluster_fasta_path, log_file_path)
    freq_df, reads_in_cluster = count_base_frequencies(cluster_fasta_path, log_file_path)

    if freq_df is None or reads_in_cluster == 0:
        return None
    
    stats_filename = os.path.join(stats_dir, f"{cluster_id}_base_freqs.csv")
    try:
        os.makedirs(stats_dir, exist_ok=True)
        freq_df.to_csv(stats_filename, index=False)
    except Exception as e:
        log(f"Cluster '{cluster_id}': ERROR saving frequency stats: {e}", log_file_path)
    
    consensus_seq = find_consensus_sequence(freq_df, cluster_id, log_file_path)
    plot_data_dict = {'df': freq_df, 'id': cluster_id, 'reads': reads_in_cluster}
    return {
        'id': cluster_id, 'reads_in_entity': reads_in_cluster,
        'consensus_seq': consensus_seq, 'plot_data': plot_data_dict
    }

def analyze_single_sample_fasta(sample_fasta_path, sample_id, stats_dir, log_file_path):
    freq_df, reads_in_sample = count_base_frequencies(sample_fasta_path, log_file_path)

    if freq_df is None or reads_in_sample < MIN_READS_PER_SAMPLE_NON_POOL:
        return None

    stats_filename = os.path.join(stats_dir, f"{sample_id}_base_freqs.csv")
    try:
        os.makedirs(stats_dir, exist_ok=True)
        freq_df.to_csv(stats_filename, index=False)
    except Exception as e:
        log(f"Sample '{sample_id}': ERROR saving frequency stats: {e}", log_file_path)

    consensus_seq = find_consensus_sequence(freq_df, sample_id, log_file_path)
    plot_data_dict = {'df': freq_df, 'id': sample_id, 'reads': reads_in_sample}
    return {
        'id': sample_id,
        'reads_in_entity': reads_in_sample,
        'consensus_seq': consensus_seq,
        'plot_data': plot_data_dict
    }

def analyze_single_sample_wrapper(path_id_tuple, stats_dir_param, log_file_path_param):
    fasta_path, s_id = path_id_tuple
    try:
        return analyze_single_sample_fasta(fasta_path, s_id, stats_dir_param, log_file_path_param)
    except Exception as e:
        log(f"Error analyzing sample '{s_id}' in child process: {e}", log_file_path_param)
        return None

# ====================== Pipeline Step Functions ======================

def setup_output_directories(output_dir, subdirs, log_file_path):
    log("Setting up output directories...", log_file_path)
    try:
        os.makedirs(output_dir, exist_ok=True)
        for folder_path in subdirs:
             os.makedirs(folder_path, exist_ok=True)
        return True
    except OSError as e:
        log(f"CRITICAL Error creating output directories: {e}. Exiting.", log_file_path)
        return False

def find_input_files(input_dir, sample_suffix, log_file_path):
    if not os.path.isdir(input_dir):
        log(f"Error: Input dir '{input_dir}' not found. Exiting.", log_file_path); return None
    try:
        input_files_list = sorted([ 
            os.path.join(input_dir, f) for f in os.listdir(input_dir) 
            if f.endswith(sample_suffix) and os.path.isfile(os.path.join(input_dir, f)) 
        ])
    except Exception as e:
        log(f"Error: Error scanning input directory '{input_dir}': {e}. Exiting.", log_file_path); return None
    
    if not input_files_list:
        log(f"Error: No files found matching '{sample_suffix}' in '{input_dir}'. Exiting.", log_file_path); return None
    log(f"Input Scan: Found {len(input_files_list)} input file(s).", log_file_path)
    return input_files_list

def run_pooling_and_filtering_step(input_files, filtered_pool_file, seq_range, max_ee, log_file_path):
    start_time = time.time()
    quoted_input_files = [f"'{p}'" for p in input_files]
    use_zcat = input_files[0].endswith(".gz") if input_files else False
    cat_cmd = "zcat -f" if use_zcat else "cat"
    
    filter_cmd = ( f"{cat_cmd} {' '.join(quoted_input_files)} | "
                   f"seqkit subseq --quiet -r {seq_range} | "
                   f"vsearch --fastq_filter - --fastq_maxee {max_ee} --fastaout '{filtered_pool_file}' --quiet" )
    
    filter_result = run_command(filter_cmd, log_file_path)
    if filter_result is None or filter_result.returncode != 0:
        log(f"Error: Pooled read filtering pipeline failed. Exiting.", log_file_path); return None, 0
    
    if not os.path.exists(filtered_pool_file) or os.path.getsize(filtered_pool_file) == 0:
        log(f"Error: Pooled filtering produced empty/missing output: {filtered_pool_file}. Exiting.", log_file_path); return None, 0
        
    pooled_read_count = count_sequences_fasta(filtered_pool_file, log_file_path)
    log(f"Pooling/Filtering: Complete ({time.time() - start_time:.1f}s). Kept {pooled_read_count:,} reads.", log_file_path)
    return filtered_pool_file, pooled_read_count

def run_clustering_step(filtered_pool_file, cluster_id_threshold, centroids_file, cluster_file_prefix, num_threads, log_file_path):
    start_time = time.time()
    clusters_dir = os.path.dirname(cluster_file_prefix)
    os.makedirs(clusters_dir, exist_ok=True)

    cluster_cmd = ( f"vsearch --cluster_fast '{filtered_pool_file}' --id {cluster_id_threshold} "
                    f"--centroids '{centroids_file}' --clusters '{cluster_file_prefix}' "
                    f"--threads {num_threads} --quiet" )
    
    cluster_result = run_command(cluster_cmd, log_file_path)
    if cluster_result is None or cluster_result.returncode != 0:
        log(f"Error: VSEARCH clustering failed. Exiting.", log_file_path); return None

    base_prefix = os.path.basename(cluster_file_prefix)
    glob_pattern = os.path.join(clusters_dir, f"{base_prefix}*")
    initial_cluster_files = [
        f for f in glob.glob(glob_pattern) 
        if os.path.basename(f).replace(base_prefix, "").replace(".fasta","").replace(".fa","").isdigit()
    ]
    initial_cluster_files.sort()

    if not initial_cluster_files:
        log(f"Error: Clustering command run, but no cluster files (e.g., {base_prefix}0) found. Exiting.", log_file_path); return None
    
    log(f"Clustering: Complete ({time.time() - start_time:.1f}s). Found {len(initial_cluster_files)} initial cluster files.", log_file_path)
    return initial_cluster_files

def filter_clusters_by_size(initial_cluster_files, min_reads_cluster_threshold, log_file_path):
    valid_cluster_files = []
    removed_count = 0
    log(f"Cluster Size Filtering: Processing {len(initial_cluster_files)} initial clusters (min reads: {min_reads_cluster_threshold})...", log_file_path)
    start_time = time.time()

    for i, cluster_file_path in enumerate(initial_cluster_files):
        reads_in_cluster = count_sequences_fasta(cluster_file_path, log_file_path)
        if reads_in_cluster >= min_reads_cluster_threshold:
            valid_cluster_files.append(cluster_file_path)
        else:
            try: os.remove(cluster_file_path)
            except OSError: pass
            removed_count += 1
        if (i + 1) % LOG_INCREMENT_LOOPS == 0 or (i + 1) == len(initial_cluster_files):
             log(f"Cluster Size Filtering: Processed {i+1}/{len(initial_cluster_files)} clusters...", log_file_path, console_output=False)
    
    log(f"Cluster Size Filtering: Complete ({time.time() - start_time:.1f}s). Removed {removed_count}, Kept {len(valid_cluster_files)} clusters.", log_file_path)
    if not valid_cluster_files:
        log(f"Error: No clusters remaining after size filtering. Exiting.", log_file_path); return None
    return valid_cluster_files

def filter_single_sample_reads(sample_fastq_path, sample_id, filtered_samples_dir, seq_range, max_ee, log_file_path):
    os.makedirs(filtered_samples_dir, exist_ok=True)
    output_fasta = os.path.join(filtered_samples_dir, f"{sample_id}_filtered.fasta")

    use_zcat = sample_fastq_path.endswith(".gz")
    cat_cmd = "zcat -f" if use_zcat else "cat"
    filter_cmd = (
        f"{cat_cmd} '{sample_fastq_path}' | "
        f"seqkit subseq --quiet -r {seq_range} | "
        f"vsearch --fastq_filter - --fastq_maxee {max_ee} --fastaout '{output_fasta}' --quiet"
    )
    filter_result = run_command(filter_cmd, log_file_path)

    if filter_result is None or filter_result.returncode != 0:
        log(f"Error: Filtering command failed for sample '{sample_id}'.", log_file_path)
        return None, 0

    read_count = count_sequences_fasta(output_fasta, log_file_path)
    if read_count == 0:
        try: os.remove(output_fasta)
        except OSError: pass
        return None, 0
    return output_fasta, read_count

def analyze_entities_parallel(iterable_of_args, process_function_partial, entity_type_desc, num_workers, log_file_path):
    num_entities_to_analyze = len(iterable_of_args)
    log(f"Parallel Analysis: Submitting {num_entities_to_analyze} {entity_type_desc.lower()} tasks ({num_workers} workers)...", log_file_path)
    
    analysis_results_dict = {}
    futures_list = []
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for args_item in iterable_of_args:
            future = executor.submit(process_function_partial, args_item)
            futures_list.append(future)

        processed_count = 0
        failed_count = 0
        for future in as_completed(futures_list):
            try:
                result = future.result()
                if result is not None:
                    analysis_results_dict[future] = result
                else:
                    failed_count += 1
            except Exception as e:
                log(f"Error (Parallel Analysis): Task failed with exception: {e}", log_file_path)
                failed_count += 1
            finally:
                processed_count += 1
                if processed_count % LOG_INCREMENT_PARALLEL == 0 or processed_count == num_entities_to_analyze:
                    log(f"Parallel Analysis: Processed {processed_count}/{num_entities_to_analyze} {entity_type_desc.lower()}...", log_file_path, console_output=False)
    
    analysis_results_list = [analysis_results_dict[f] for f in futures_list if f in analysis_results_dict]
    
    log(f"Parallel Analysis: {entity_type_desc} analysis finished. Successful: {len(analysis_results_list)}, Failed/Skipped: {failed_count + (num_entities_to_analyze - processed_count)}.", log_file_path)

    if not analysis_results_list and num_entities_to_analyze > 0 :
        log(f"Error: No {entity_type_desc.lower()} successfully analyzed. Exiting.", log_file_path); return None, []

    entity_consensus_results = [(res['id'], res['consensus_seq']) for res in analysis_results_list if res.get('consensus_seq')]
    log(f"Parallel Analysis: Generated {len(entity_consensus_results)} non-empty consensus sequences.", log_file_path)
    return analysis_results_list, entity_consensus_results

def generate_plots_for_entities(analysis_results_list, plots_dir, log_file_path, entity_type_desc):
    if not analysis_results_list: return 0,0
    
    num_results = len(analysis_results_list)
    log(f"Plot Generation: Generating plots for {num_results} {entity_type_desc.lower()} (Serial)...", log_file_path)
    plots_generated_count = 0
    plots_failed_count = 0
    
    for i, result_data in enumerate(analysis_results_list):
        plot_info = result_data.get('plot_data')
        if (plot_info and isinstance(plot_info.get('df'), pd.DataFrame) and
            not plot_info['df'].empty and 'id' in plot_info and 'reads' in plot_info):
            
            if plot_base_frequencies(plot_info['df'], plot_info['id'], plot_info['reads'], plots_dir, log_file_path):
                plots_generated_count += 1
            else: plots_failed_count += 1
        else:
            plots_failed_count += 1
        
        if (i + 1) % LOG_INCREMENT_LOOPS == 0 or (i + 1) == num_results:
             log(f"Plot Generation: Processed {i+1}/{num_results} plots...", log_file_path, console_output=False)

    log(f"Plot Generation: Finished. Generated {plots_generated_count}, Failed/Skipped {plots_failed_count} for {entity_type_desc.lower()}.", log_file_path)
    return plots_generated_count, plots_failed_count

def run_blast_pipeline(entity_consensus_results, blast_input_fasta_name, blast_db_path, blast_evalue_param, blast_outfmt_param, blast_output_tsv_path, best_hits_csv_path, num_threads_param, log_file_path, entity_type_desc, analysis_results_list=None):
    if not entity_consensus_results:
        log(f"BLAST: Skipping for {entity_type_desc} - No consensus sequences.", log_file_path)
        return False

    log(f"BLAST: Preparing and Running BLAST for {len(entity_consensus_results)} {entity_type_desc.lower()} consensus sequences...", log_file_path)

    blast_dir = os.path.dirname(blast_output_tsv_path)
    os.makedirs(blast_dir, exist_ok=True)
    blast_input_fasta_path = os.path.join(blast_dir, blast_input_fasta_name)

    try:
        with open(blast_input_fasta_path, "w") as f_out:
            for entity_id, sequence in entity_consensus_results:
                f_out.write(f">{entity_id}\n{sequence}\n")
    except IOError as e:
        log(f"Error (BLAST): Error writing BLAST input {blast_input_fasta_path}: {e}. Skipping BLAST.", log_file_path)
        return False

    log(f"BLAST: Running BLASTN vs '{blast_db_path}' (e-val<{blast_evalue_param})...", log_file_path)
    max_len = max(len(s) for _, s in entity_consensus_results) if entity_consensus_results else 0
    task = "blastn-short" if 0 < max_len < 50 else "blastn"
    cmd = (f"blastn -query '{blast_input_fasta_path}' -db '{blast_db_path}' -out '{blast_output_tsv_path}' "
           f"-outfmt \"{blast_outfmt_param}\" -evalue {blast_evalue_param} -task {task} -num_threads {num_threads_param} -strand plus ")
    
    result = run_command(cmd, log_file_path)
    if result is None or result.returncode != 0:
        log(f"Warning (BLAST): BLAST command may have failed. Output will be processed if available.", log_file_path)

    log(f"BLAST: Processing BLAST output for {entity_type_desc}...", log_file_path)
    processed_successfully = False

    all_submitted_queries_data = []
    entity_id_to_reads_map = {}
    if analysis_results_list:
        entity_id_to_reads_map = {
            item['id']: item['reads_in_entity'] 
            for item in analysis_results_list if 'id' in item and 'reads_in_entity' in item
        }

    for entity_id, _ in entity_consensus_results:
        reads_count = entity_id_to_reads_map.get(entity_id, 0) 
        all_submitted_queries_data.append({'query_id': entity_id, 'query_reads_count': reads_count})
    
    if not all_submitted_queries_data:
        log(f"BLAST: No query data to process for {entity_type_desc.lower()}.", log_file_path)
        return False
            
    all_queries_df = pd.DataFrame(all_submitted_queries_data)

    hits_from_blast_df = pd.DataFrame() 
    blast_field_names_from_outfmt = blast_outfmt_param.split(" ")[1:] 
    csv_column_names = [BLAST_HEADER_MAP.get(h, h) for h in blast_field_names_from_outfmt]

    if os.path.exists(blast_output_tsv_path) and os.path.getsize(blast_output_tsv_path) > 0:
        try:
            df_raw_blast = pd.read_csv(blast_output_tsv_path, sep="\t", names=csv_column_names)
            if not df_raw_blast.empty:
                if 'bit_score' in df_raw_blast.columns and 'query_id' in df_raw_blast.columns:
                    df_raw_blast['bit_score'] = pd.to_numeric(df_raw_blast['bit_score'], errors='coerce')
                    df_raw_blast.dropna(subset=['bit_score'], inplace=True)
                    if not df_raw_blast.empty:
                        hits_from_blast_df = df_raw_blast.sort_values(
                            ["query_id", "bit_score"], ascending=[True, False]
                        )
                else:
                    missing_cols = [col for col in ['bit_score', 'query_id'] if col not in df_raw_blast.columns]
                    log(f"Warning (BLAST): Key column(s) {missing_cols} not found in raw BLAST output. Cannot process hits.", log_file_path)
        except Exception as e:
            log(f"Error (BLAST): Failed to process raw BLAST output {os.path.basename(blast_output_tsv_path)}: {e}", log_file_path)

    if not all_queries_df.empty:
        if not hits_from_blast_df.empty:
            final_df = pd.merge(all_queries_df, hits_from_blast_df, on="query_id", how="left")
        else: 
            final_df = all_queries_df.copy()
            for col_name in csv_column_names:
                if col_name not in final_df.columns:
                    final_df[col_name] = pd.NA 
        
        for col in final_df.columns:
            final_df[col] = final_df[col].astype("object").fillna("N/A").infer_objects(copy=False)

        desired_column_order = ['query_id', 'query_reads_count']
        for mapped_col in csv_column_names:
            if mapped_col != 'query_id' and mapped_col not in desired_column_order:
                desired_column_order.append(mapped_col)
        
        for col in desired_column_order:
            if col not in final_df.columns:
                final_df[col] = "N/A"
        final_df = final_df[desired_column_order]

        try:
            final_df.to_csv(best_hits_csv_path, index=False)
            log(f"BLAST: Final report for {len(final_df)} query-hit pairs saved: {os.path.basename(best_hits_csv_path)}", log_file_path)
            processed_successfully = True
        except Exception as e:
            log(f"Error (BLAST): Failed to save final BLAST report to {os.path.basename(best_hits_csv_path)}: {e}", log_file_path)
    else: 
        log(f"BLAST: No submitted queries to generate a final report for {entity_type_desc}.", log_file_path)

    return processed_successfully

# ====================== Main Execution ======================

def main():
    global THREADS, BLAST_EVALUE 

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", dest="input_dir", required=True, help="Input folder with FASTQ files.")
    parser.add_argument("-o", "--output", dest="output_dir", required=True, help="Output folder for results.")
    parser.add_argument("-s", "--suffix", dest="sample_suffix", required=True, help="File suffix of target files (e.g., _R1.fastq.gz).")
    parser.add_argument("-p", "--pool", action="store_true", help="Pool samples for joint filtering/clustering (default: per-sample).")
    parser.add_argument("-t", "--threads", type=int, default=THREADS, help="Number of threads/processes.")
    parser.add_argument("-e", "--e_value", type=float, default=BLAST_EVALUE, help="E-value threshold for BLAST.")
    args = parser.parse_args()

    input_dir, output_dir, sample_suffix = args.input_dir, args.output_dir, args.sample_suffix
    THREADS = args.threads
    BLAST_EVALUE = args.e_value
    run_pooled_mode = args.pool

    pipeline_start_time = time.time()
    logs_dir = os.path.join(output_dir, "logs")
    log_file = os.path.join(logs_dir, "pipeline.log")
    if not setup_logging(log_file): sys.exit(1)

    blast_db_to_use = os.path.join(os.environ.get('APHUNTER_BLASTDB_DIR', '.'), 
                                   os.environ.get('APHUNTER_BLASTDB_NAME', 'its_primer_db'))
    db_source_msg = f"'{blast_db_to_use}' (from APHUNTER_BLASTDB_DIR/NAME env vars if set, else defaults)"
    if not (os.environ.get('APHUNTER_BLASTDB_DIR') and os.environ.get('APHUNTER_BLASTDB_NAME')):
        log(f"Warning: APHUNTER_BLASTDB_DIR or NAME env vars not set. Using default/relative path for DB: {blast_db_to_use}", log_file)
    
    log("================ Pipeline Start ================", log_file)
    log(f"Input Dir: {os.path.abspath(input_dir)}", log_file)
    log(f"Output Dir: {os.path.abspath(output_dir)}", log_file)
    log(f"Run Mode: {'Pooled Clustering' if run_pooled_mode else 'Per-Sample Analysis'}", log_file)
    log(f"Threads: {THREADS}", log_file)
    log(f"BLAST DB: {db_source_msg}", log_file, console_output=False)
    log(f"BLAST E-value: {BLAST_EVALUE}", log_file, console_output=False)

    plots_dir = os.path.join(output_dir, "plots")
    stats_dir = os.path.join(output_dir, "stats")
    blast_dir = os.path.join(output_dir, "blast")
    
    analysis_results_list = []
    entity_consensus_results = []
    final_blast_output_csv = ""
    initial_sample_count = 0

    common_subdirs = [plots_dir, stats_dir, blast_dir, logs_dir]
    if run_pooled_mode:
        clusters_dir = os.path.join(output_dir, "clusters")
        all_subdirs = common_subdirs + [clusters_dir]
    else:
        filtered_samples_output_dir = os.path.join(output_dir, "filtered_samples")
        all_subdirs = common_subdirs + [filtered_samples_output_dir]
    if not setup_output_directories(output_dir, all_subdirs, log_file): sys.exit(1)

    step_num = 1
    log(f"Step {step_num}: Scanning input files...", log_file)
    input_files = find_input_files(input_dir, sample_suffix, log_file)
    if input_files is None: sys.exit(1)
    initial_sample_count = len(input_files)
    step_num += 1
    
    final_entities_count = 0
    final_entities_analyzed_count = 0
    final_entities_for_blast_count = 0

    if run_pooled_mode:
        log(f"--- Starting POOLED MODE Pipeline ---", log_file)
        filtered_pool_file = os.path.join(output_dir, "pooled_filtered.fasta")
        all_centroids_file = os.path.join(output_dir, "all_centroids.fasta")
        cluster_file_prefix = os.path.join(clusters_dir, "cluster.")

        log(f"Step {step_num}: Pooling & Filtering reads...", log_file)
        filtered_pool_file, pooled_read_count = run_pooling_and_filtering_step(input_files, filtered_pool_file, SEQKIT_RANGE, MAX_EE, log_file)
        if filtered_pool_file is None: sys.exit(1)
        step_num += 1
        
        min_reads_cluster = max(MIN_READS_FLOOR, pooled_read_count // 2000 if pooled_read_count > 0 else 1)

        log(f"Step {step_num}: Clustering reads (ID >= {CLUSTER_IDENTITY_THRESHOLD})...", log_file)
        initial_cluster_files = run_clustering_step(filtered_pool_file, CLUSTER_IDENTITY_THRESHOLD, all_centroids_file, cluster_file_prefix, THREADS, log_file)
        if initial_cluster_files is None: sys.exit(1)
        try: os.remove(filtered_pool_file)
        except OSError: pass
        step_num += 1

        log(f"Step {step_num}: Filtering clusters by size (Min reads: {min_reads_cluster})...", log_file)
        valid_cluster_files = filter_clusters_by_size(initial_cluster_files, min_reads_cluster, log_file)
        if valid_cluster_files is None: sys.exit(1)
        final_entities_count = len(valid_cluster_files)
        step_num += 1

        if final_entities_count > 0:
            log(f"Step {step_num}: Analyzing {final_entities_count} clusters...", log_file)
            process_func_for_pool = functools.partial(process_cluster_file, stats_dir=stats_dir, log_file_path=log_file)
            analysis_results_list, entity_consensus_results = analyze_entities_parallel(valid_cluster_files, process_func_for_pool, "Clusters", THREADS, log_file)
            if analysis_results_list is None: sys.exit(1)
            final_entities_analyzed_count = len(analysis_results_list)
            final_entities_for_blast_count = len(entity_consensus_results)
            step_num += 1

            log(f"Step {step_num}: Generating plots for {final_entities_analyzed_count} clusters...", log_file)
            generate_plots_for_entities(analysis_results_list, plots_dir, log_file, "Clusters")
            step_num += 1

            log(f"Step {step_num}: Running BLAST for {final_entities_for_blast_count} cluster consensus...", log_file)
            blast_output_tsv_pooled = os.path.join(blast_dir, "pooled_blast_output.tsv")
            final_blast_output_csv = os.path.join(blast_dir, "pooled_final_blast_results.csv")
            run_blast_pipeline(entity_consensus_results, "cluster_consensus_for_blast.fasta", blast_db_to_use, BLAST_EVALUE, BLAST_OUTFMT, blast_output_tsv_pooled, final_blast_output_csv, THREADS, log_file, "Clusters", analysis_results_list)
    
    else: # NON-POOLED (PER-SAMPLE) MODE
        log(f"--- Starting NON-POOLED (Per-Sample) MODE Pipeline ---", log_file)
        
        log(f"Step {step_num}: Filtering individual samples (Min reads: {MIN_READS_PER_SAMPLE_NON_POOL})...", log_file)
        filtered_sample_fastas_args = []
        samples_discarded_count = 0
        for i, fq_path in enumerate(input_files):
            sample_id = get_sample_id_from_path(fq_path, sample_suffix)
            filtered_fasta, read_count = filter_single_sample_reads(fq_path, sample_id, filtered_samples_output_dir, SEQKIT_RANGE, MAX_EE, log_file)
            if filtered_fasta and read_count >= MIN_READS_PER_SAMPLE_NON_POOL:
                filtered_sample_fastas_args.append((filtered_fasta, sample_id))
            else: samples_discarded_count += 1
            
            if (i + 1) % LOG_INCREMENT_LOOPS == 0 or (i + 1) == len(input_files):
                 log(f"Sample Filtering: Processed {i+1}/{len(input_files)} samples...", log_file, console_output=False)

        final_entities_count = len(filtered_sample_fastas_args)
        log(f"Sample Filtering: Complete. Kept {final_entities_count}, discarded {samples_discarded_count} samples.", log_file)
        if not filtered_sample_fastas_args:
            log(f"Error: No samples remaining after filtering. Exiting.", log_file); sys.exit(1)
        step_num += 1

        log(f"Step {step_num}: Analyzing {final_entities_count} filtered samples...", log_file)
        analyze_sample_func_partial = functools.partial(analyze_single_sample_wrapper,
                                                        stats_dir_param=stats_dir,
                                                        log_file_path_param=log_file)
        analysis_results_list, entity_consensus_results = analyze_entities_parallel(
            filtered_sample_fastas_args, analyze_sample_func_partial, "Samples", THREADS, log_file
        )
        if analysis_results_list is None: sys.exit(1)
        final_entities_analyzed_count = len(analysis_results_list)
        final_entities_for_blast_count = len(entity_consensus_results)
        step_num += 1

        log(f"Step {step_num}: Generating plots for {final_entities_analyzed_count} samples...", log_file)
        generate_plots_for_entities(analysis_results_list, plots_dir, log_file, "Samples")
        step_num += 1

        log(f"Step {step_num}: Running BLAST for {final_entities_for_blast_count} sample consensus...", log_file)
        blast_output_tsv_nonpool = os.path.join(blast_dir, "per_sample_blast_output.tsv")
        final_blast_output_csv = os.path.join(blast_dir, "per_sample_final_blast_results.csv")
        run_blast_pipeline(entity_consensus_results, "sample_consensus_for_blast.fasta", blast_db_to_use, BLAST_EVALUE, BLAST_OUTFMT, blast_output_tsv_nonpool, final_blast_output_csv, THREADS, log_file, "Samples", analysis_results_list)

    elapsed_seconds = time.time() - pipeline_start_time
    elapsed_minutes, elapsed_sec_rem = divmod(elapsed_seconds, 60)
    log(f"==================== Pipeline Summary ====================", log_file)
    log(f"Pipeline Finished. Total Time: {int(elapsed_minutes)}m {elapsed_sec_rem:.1f}s.", log_file)
    log(f"Output Directory: {os.path.abspath(output_dir)}", log_file)
    log(f"Initial Input Samples: {initial_sample_count}", log_file)

    entity_label = "Clusters" if run_pooled_mode else "Samples"
    log(f"{entity_label} Kept for Analysis: {final_entities_count}", log_file)
    log(f"{entity_label} Successfully Analyzed: {final_entities_analyzed_count}", log_file)
    log(f"{entity_label} Consensus for BLAST: {final_entities_for_blast_count}", log_file)

    if entity_consensus_results and os.path.exists(final_blast_output_csv) and os.path.getsize(final_blast_output_csv) > 0:
        try:
            with open(final_blast_output_csv, 'r') as f_blast_csv:
                header_line = f_blast_csv.readline().strip()
                actual_headers_in_file = header_line.split(',')
                if 'query_reads_count' in actual_headers_in_file:
                    log(f"Info: 'query_reads_count' column confirmed in {os.path.basename(final_blast_output_csv)}.", log_file, console_output=False)
                else:
                    log(f"Warning: 'query_reads_count' column NOT found in {os.path.basename(final_blast_output_csv)} header.", log_file)
            
            with open(final_blast_output_csv, 'r') as f_blast_csv_count:
                num_total_entries = max(0, sum(1 for _ in f_blast_csv_count) -1) # -1 for header
            log(f"Final BLAST Report Entries (query-hit pairs): {num_total_entries} (in {os.path.basename(final_blast_output_csv)})", log_file)
        except Exception as e_count:
            log(f"Final BLAST Report: {os.path.basename(final_blast_output_csv)} (exists, but count/header check failed: {e_count})", log_file)
    elif entity_consensus_results:
        log(f"BLAST was run, but final report '{os.path.basename(final_blast_output_csv)}' is missing or empty.", log_file)
    else:
        log(f"BLAST step was skipped (no consensus sequences).", log_file)
    
    log(f"Detailed logs: {log_file}", log_file)
    log(f"==========================================================", log_file)

if __name__ == "__main__":
    main()