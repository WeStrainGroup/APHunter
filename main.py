#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import collections
import atexit
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
import csv
import glob
import functools
import shlex
import sys
import threading

# ====================== Parameter Settings ======================
THREADS = max(8, multiprocessing.cpu_count() // 2 if multiprocessing.cpu_count() else 16)
SEQKIT_RANGE = "1:100"
MIN_READS_FLOOR = 100
MAX_EE = 1.0
MIN_READS_PER_SAMPLE_NON_POOL = 1000
LOG_INCREMENT_LOOPS = 25
LOG_INCREMENT_PARALLEL = 50
CLUSTER_IDENTITY_THRESHOLD = 0.9
# Step-prefixed intermediate folders; final BLAST CSVs go under OUT_FINAL only.
OUT_NONPOOL_FILTERED = "1_filtered_samples"
OUT_NONPOOL_STATS = "2_stats"
OUT_NONPOOL_PLOTS = "3_plots"
OUT_NONPOOL_BLAST_WORK = "4_blast_work"
OUT_POOL_WORK = "1_pool_work"
OUT_POOL_CLUSTERS = "2_clusters"
OUT_POOL_STATS = "3_stats"
OUT_POOL_PLOTS = "4_plots"
OUT_POOL_BLAST_WORK = "5_blast_work"
OUT_FINAL = "final"

BLAST_EVALUE = 1e-3
BLAST_WORD_SIZE = 7
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

# Compact, stable CSV headers for BLAST summary files (query/hit first, provenance last).
BLAST_EXPORT_RENAME = {
    "blast_db_prefix": "db_path",
    "blast_db_label": "db",
    "blast_db_resolution": "db_from",
    "query_id": "query",
    "query_reads_count": "reads",
    "hit_rank": "rank",
    "top_suspect_recommendation": "best_hit",
    "%_identity": "pct_id",
    "alignment_length": "aln_len",
    "subject_id": "subject",
    "mismatches": "mismatch",
    "gap_opens": "gaps",
    "query_length": "qlen",
    "q_start": "qstart",
    "q_end": "qend",
    "subject_length": "slen",
    "s_start": "sstart",
    "s_end": "send",
}
BLAST_EXPORT_LEADING = (
    "query",
    "reads",
    "rank",
    "best_hit",
    "subject",
    "pct_id",
    "aln_len",
    "evalue",
    "bit_score",
)
BLAST_EXPORT_TRAILING = ("db", "db_from", "db_path")


def prepare_blast_results_for_export(df):
    """Rename to short headers and reorder columns for scanning."""
    renamed = df.rename(columns={k: v for k, v in BLAST_EXPORT_RENAME.items() if k in df.columns})
    cols = list(renamed.columns)
    leading = [c for c in BLAST_EXPORT_LEADING if c in cols]
    trailing = [c for c in BLAST_EXPORT_TRAILING if c in cols]
    middle = [c for c in cols if c not in leading and c not in trailing]
    return renamed[leading + middle + trailing]


TOP_SUSPECT_EXPORT_COLS = (
    "query",
    "reads",
    "status",
    "subject",
    "pct_id",
    "aln_len",
    "evalue",
    "bit_score",
    "db",
    "db_from",
)

PLOT_DPI = 300
COLOR_MAP = {'A': '#FF6347', 'T': '#4682B4', 'C': '#008000', 'G': '#800080'}
THRESHOLD_FREQ = 0.75
CONSENSUS_WARNING_CONSOLE = False
IUPAC_BASE_WEIGHTS = {
    'A': {'A': 1.0},
    'T': {'T': 1.0},
    'C': {'C': 1.0},
    'G': {'G': 1.0},
    # IUPAC ambiguity codes: split a single read vote across candidates.
    'R': {'A': 0.5, 'G': 0.5},
    'Y': {'C': 0.5, 'T': 0.5},
    'S': {'G': 0.5, 'C': 0.5},
    'W': {'A': 0.5, 'T': 0.5},
    'K': {'G': 0.5, 'T': 0.5},
    'M': {'A': 0.5, 'C': 0.5},
    'B': {'C': 1.0 / 3.0, 'G': 1.0 / 3.0, 'T': 1.0 / 3.0},
    'D': {'A': 1.0 / 3.0, 'G': 1.0 / 3.0, 'T': 1.0 / 3.0},
    'H': {'A': 1.0 / 3.0, 'C': 1.0 / 3.0, 'T': 1.0 / 3.0},
    'V': {'A': 1.0 / 3.0, 'C': 1.0 / 3.0, 'G': 1.0 / 3.0},
}

# BLAST nucleotide DB: any of these files adjacent to the prefix indicates a valid DB.
_BLAST_NT_INDEX_SUFFIXES = (
    '.nhr', '.nin', '.nsq',
    '.ndb', '.not', '.ntf', '.nto', '.nos',
    '.njs',
    '.nhd', '.nti', '.nnd', '.nni',
    '.nal',
)
# Centralized logging queue/listener for multi-process safe writes.
LOG_QUEUE = None
LOG_LISTENER_THREAD = None
LOG_FILE_PATH_GLOBAL = None


def _project_root_dir():
    return os.path.dirname(os.path.abspath(__file__))


def _blast_cli_builtin_db(s):
    """Normalize --blast-db choice (case-insensitive)."""
    return str(s).strip().lower()


def blast_nucleotide_db_prefix_exists(prefix):
    p = os.path.abspath(os.path.expanduser(prefix))
    return any(os.path.isfile(p + ext) for ext in _BLAST_NT_INDEX_SUFFIXES)


def resolve_blast_database(args):
    """
    Resolve BLAST nucleotide DB prefix for `blastn -db`.

    Priority:
      1) --blast-db-dir (custom directory; auto-detect one prefix; label = basename of prefix)
      2) --blast-db (built-in: its | 16s). Mutually exclusive with --blast-db-dir at CLI parse time.
      3) ENV fallback: APHUNTER_BLASTDB_DIR + APHUNTER_BLASTDB_NAME
      4) default built-in: its
    Returns (prefix_or_None, meta_dict). On failure meta_dict contains 'error'.
    """
    built_in_db_map = {
        'its': os.path.join(_project_root_dir(), 'blast_db', 'its_primer_db', 'its_primers'),
        '16s': os.path.join(_project_root_dir(), 'blast_db', '16s_primer_db', '16s_primers'),
    }

    requested_dir = (os.path.abspath(os.path.expanduser(str(getattr(args, 'blast_db_dir')).strip()))
                     if getattr(args, 'blast_db_dir', None) and str(getattr(args, 'blast_db_dir')).strip()
                     else None)

    if requested_dir:
        env_name = os.environ.get('APHUNTER_BLASTDB_NAME')
        if env_name and str(env_name).strip():
            env_prefix = os.path.join(requested_dir, str(env_name).strip())
            if blast_nucleotide_db_prefix_exists(env_prefix):
                return env_prefix, {
                    'prefix': env_prefix,
                    'label': os.path.basename(env_prefix),
                    'resolution': 'CLI --blast-db-dir + ENV APHUNTER_BLASTDB_NAME',
                }

        matched_prefixes = set()
        for ext in _BLAST_NT_INDEX_SUFFIXES:
            for file_path in glob.glob(os.path.join(requested_dir, f"*{ext}")):
                matched_prefixes.add(file_path[:-len(ext)])
        matched_prefixes = sorted(matched_prefixes)

        if len(matched_prefixes) == 1:
            prefix = matched_prefixes[0]
            return prefix, {
                'prefix': prefix,
                'label': os.path.basename(prefix),
                'resolution': 'CLI --blast-db-dir (auto-detected prefix)',
            }

        if len(matched_prefixes) > 1:
            prefix_names = ", ".join(os.path.basename(p) for p in matched_prefixes[:6])
            suffix_msg = " ..." if len(matched_prefixes) > 6 else ""
            return None, {
                'error': (
                    f"Multiple BLAST DB prefixes found in {requested_dir!r}: {prefix_names}{suffix_msg}. "
                    f"Please keep only one BLAST nucleotide index set in this directory."
                ),
                'prefix': requested_dir,
                'label': os.path.basename(requested_dir),
                'resolution': 'CLI --blast-db-dir',
            }

        return None, {
            'error': f"No BLAST nucleotide index files found in directory {requested_dir!r}.",
            'prefix': requested_dir,
            'label': os.path.basename(requested_dir),
            'resolution': 'CLI --blast-db-dir',
        }

    requested_builtin = getattr(args, 'blast_db', None)
    if requested_builtin:
        name = str(requested_builtin).strip().lower()
        if name not in built_in_db_map:
            return None, {
                'error': f"Unsupported --blast-db value: {name!r}. Use 'its' or '16s', or use --blast-db-dir for a custom index directory.",
                'prefix': '',
                'label': name,
                'resolution': 'CLI --blast-db',
            }
        prefix = built_in_db_map[name]
        if not blast_nucleotide_db_prefix_exists(prefix):
            return None, {
                'error': f"Built-in database '{name}' not found at prefix {prefix!r}.",
                'prefix': prefix,
                'label': name,
                'resolution': 'CLI --blast-db',
            }
        return prefix, {
            'prefix': prefix,
            'label': name,
            'resolution': 'CLI --blast-db',
        }

    env_dir = os.environ.get('APHUNTER_BLASTDB_DIR')
    env_name = os.environ.get('APHUNTER_BLASTDB_NAME')
    if env_dir and env_name:
        prefix = os.path.join(os.path.abspath(os.path.expanduser(env_dir.strip())), env_name.strip())
        if blast_nucleotide_db_prefix_exists(prefix):
            return prefix, {
                'prefix': prefix,
                'label': env_name.strip(),
                'resolution': 'ENV APHUNTER_BLASTDB_DIR + APHUNTER_BLASTDB_NAME',
            }

    prefix = built_in_db_map['its']
    if not blast_nucleotide_db_prefix_exists(prefix):
        return None, {
            'error': f"Default built-in ITS database not found at prefix {prefix!r}.",
            'prefix': prefix,
            'label': 'its',
            'resolution': 'default built-in (its)',
        }
    return prefix, {
        'prefix': prefix,
        'label': 'its',
        'resolution': 'default built-in (its)',
    }


def sanitize_db_label_for_filename(label):
    """Normalize DB label for safe output filenames."""
    raw = str(label).strip() if label is not None else ""
    if not raw:
        raw = "unknown_db"
    safe = []
    for ch in raw:
        if ch.isalnum() or ch in ("-", "_"):
            safe.append(ch)
        else:
            safe.append("_")
    normalized = "".join(safe)
    while "__" in normalized:
        normalized = normalized.replace("__", "_")
    normalized = normalized.strip("_")
    return normalized if normalized else "unknown_db"


def infer_read_orientation(sample_suffix):
    """Infer forward/reverse read orientation from input filename suffix."""
    suffix = str(sample_suffix or "").lower()
    reverse_markers = (".2.fq.gz", ".2.fastq.gz", ".2.fq", ".2.fastq", "_r2", ".r2.")
    forward_markers = (".1.fq.gz", ".1.fastq.gz", ".1.fq", ".1.fastq", "_r1", ".r1.")
    if any(marker in suffix for marker in reverse_markers) or suffix.endswith(".2.gz"):
        return "reverse"
    if any(marker in suffix for marker in forward_markers) or suffix.endswith(".1.gz"):
        return "forward"
    return None


def classify_primer_orientation(subject_id, db_label):
    """Classify a BLAST subject as forward or reverse primer."""
    name = str(subject_id or "").strip()
    if not name:
        return None
    lower = name.lower()
    db = str(db_label or "").lower()
    if db == "its" or "_fwd" in lower or "_rev" in lower:
        if lower.endswith("_rev") or lower.endswith("_rev_rc"):
            return "reverse"
        if lower.endswith("_fwd") or lower.endswith("_fwd_rc"):
            return "forward"
        if "_rev" in lower:
            return "reverse"
        if "_fwd" in lower:
            return "forward"
    if lower.endswith("r"):
        return "reverse"
    if lower.endswith("f"):
        return "forward"
    return None


def _rerank_blast_hits(ranked_hits_df):
    if ranked_hits_df.empty:
        return ranked_hits_df
    sort_cols = [c for c in ['query_id', 'bit_score', '%_identity', 'alignment_length', 'evalue'] if c in ranked_hits_df.columns]
    sort_ascending = []
    for col in sort_cols:
        if col == 'query_id':
            sort_ascending.append(True)
        elif col == 'evalue':
            sort_ascending.append(True)
        else:
            sort_ascending.append(False)
    ranked_hits_df = ranked_hits_df.sort_values(sort_cols, ascending=sort_ascending, na_position='last')
    ranked_hits_df['hit_rank'] = ranked_hits_df.groupby('query_id', sort=False).cumcount() + 1
    ranked_hits_df['is_top_suspect'] = ranked_hits_df['hit_rank'] == 1
    return ranked_hits_df


def filter_ranked_hits_by_orientation(ranked_hits_df, read_orientation, db_label, log_file_path):
    """Keep only primer hits matching the read orientation (R1 -> forward, R2 -> reverse)."""
    if ranked_hits_df.empty or not read_orientation or 'subject_id' not in ranked_hits_df.columns:
        return ranked_hits_df
    primer_orientations = ranked_hits_df['subject_id'].map(lambda s: classify_primer_orientation(s, db_label))
    keep_mask = primer_orientations.isna() | (primer_orientations == read_orientation)
    removed = int((~keep_mask).sum())
    if removed:
        log(
            f"BLAST: orientation filter removed {removed} hit row(s) "
            f"({read_orientation} read -> keep {read_orientation} primers only)",
            log_file_path,
        )
    filtered = ranked_hits_df[keep_mask].copy()
    return _rerank_blast_hits(filtered)


def _confidence_label(top_read_pct, hit_query_pct, median_identity):
    if hit_query_pct <= 0:
        return "NONE"
    if top_read_pct >= 90.0 and median_identity >= 90.0:
        return "HIGH"
    if top_read_pct >= 50.0 or hit_query_pct >= 50.0:
        return "MEDIUM"
    return "LOW"


def build_primer_call_summary(top_suspect_df, entity_type_desc, db_label, read_orientation):
    """Build run-level primer call text lines from the top-suspect table."""
    df = top_suspect_df.copy()
    total_queries = len(df)
    reads = pd.to_numeric(df.get("reads", 0), errors="coerce").fillna(0)
    total_reads = float(reads.sum())
    hits = df[df.get("status", "") == "hit"].copy()
    hit_queries = len(hits)
    hit_reads = float(pd.to_numeric(hits.get("reads", 0), errors="coerce").fillna(0).sum()) if hit_queries else 0.0

    query_unit = "samples" if entity_type_desc.lower().startswith("sample") else "clusters"
    read_label = "forward" if read_orientation == "forward" else "reverse" if read_orientation == "reverse" else "unknown"

    lines = [
        "APHunter primer call summary",
        "===========================",
        f"Mode: {'per-sample' if query_unit == 'samples' else 'pooled (per-cluster)'}",
        f"Read orientation: {read_label}",
        f"BLAST DB: {db_label}",
        "Primer reference: https://apvisual.netlify.app/",
        "",
    ]

    if hit_queries == 0 or total_queries == 0:
        lines.extend([
            "Likely primer: none",
            f"Overall hit rate: 0/{total_queries} {query_unit} (0.0%)",
            "Confidence: NONE",
            "",
            "Per-primer breakdown: no hits",
        ])
        return lines

    hit_reads_series = pd.to_numeric(hits["reads"], errors="coerce").fillna(0)
    hit_ids = pd.to_numeric(hits.get("pct_id", pd.NA), errors="coerce")
    breakdown = (
        hits.assign(_reads=hit_reads_series)
        .groupby("subject", as_index=False)
        .agg(query_support=("query", "count"), read_support=("_reads", "sum"))
        .sort_values(["read_support", "query_support"], ascending=False)
    )
    top = breakdown.iloc[0]
    likely_subject = str(top["subject"])
    top_query_pct = (float(top["query_support"]) / total_queries) * 100.0
    top_read_pct = (float(top["read_support"]) / hit_reads) * 100.0 if hit_reads else 0.0
    likely_ids = hit_ids[hits["subject"] == likely_subject]
    median_identity = float(likely_ids.median()) if not likely_ids.empty else float("nan")

    lines.extend([
        f"Likely primer: {likely_subject}",
        f"  Query support: {int(top['query_support'])}/{total_queries} {query_unit} ({top_query_pct:.1f}%)",
        f"  Read support (among hits): {int(top['read_support']):,}/{int(hit_reads):,} ({top_read_pct:.1f}%)",
        f"  Median identity: {median_identity:.1f}%" if pd.notna(median_identity) else "  Median identity: n/a",
    ])

    if len(breakdown) > 1:
        runner = breakdown.iloc[1]
        runner_read_pct = (float(runner["read_support"]) / hit_reads) * 100.0 if hit_reads else 0.0
        runner_subject = str(runner["subject"])
        lines.append(
            f"  Runner-up: {runner_subject} "
            f"({int(runner['query_support'])} {query_unit}, {runner_read_pct:.1f}% of hit reads)"
        )
    else:
        lines.append("  Runner-up: none")

    lines.extend([
        "",
        f"Overall hit rate: {hit_queries}/{total_queries} {query_unit} ({(hit_queries / total_queries) * 100.0:.1f}%)",
        f"Confidence: {_confidence_label(top_read_pct, (hit_queries / total_queries) * 100.0, median_identity if pd.notna(median_identity) else 0.0)}",
        "",
        "Per-primer breakdown (by hit reads):",
    ])

    for _, row in breakdown.head(10).iterrows():
        read_pct = (float(row["read_support"]) / hit_reads) * 100.0 if hit_reads else 0.0
        query_pct = (float(row["query_support"]) / total_queries) * 100.0
        subject = str(row["subject"])
        lines.append(
            f"  {subject:24}  {query_unit}={int(row['query_support']):>4} ({query_pct:5.1f}%)  "
            f"reads={read_pct:5.1f}%"
        )

    return lines


def write_primer_call_summary(final_dir, top_suspect_df, entity_type_desc, db_label, read_orientation, log_file_path):
    """Write final/primer_call_summary.txt and return concise lines for terminal/log output."""
    summary_lines = build_primer_call_summary(top_suspect_df, entity_type_desc, db_label, read_orientation)
    summary_path = os.path.join(final_dir, "primer_call_summary.txt")
    try:
        os.makedirs(final_dir, exist_ok=True)
        with open(summary_path, "w", encoding="utf-8") as f_out:
            f_out.write("\n".join(summary_lines) + "\n")
        log(f"BLAST: primer call -> {os.path.basename(summary_path)}", log_file_path)
    except Exception as e:
        log(f"Warning (BLAST): Failed to save primer call summary: {e}", log_file_path, level="WARNING")
        return None

    console_lines = []
    for line in summary_lines:
        if line in ("APHunter primer call summary", "==========================="):
            continue
        if line.strip():
            console_lines.append(line)
    console_lines.append(f"Full summary: {summary_path}")
    return console_lines

# ====================== Helper Functions ======================

def _log_listener_worker():
    """Single writer thread: consumes queued log lines and appends to log file."""
    global LOG_QUEUE, LOG_FILE_PATH_GLOBAL
    while True:
        try:
            message = LOG_QUEUE.get()
        except Exception:
            continue
        if message is None:  # sentinel for shutdown
            break
        try:
            log_dir = os.path.dirname(LOG_FILE_PATH_GLOBAL)
            if log_dir:
                os.makedirs(log_dir, exist_ok=True)
            with open(LOG_FILE_PATH_GLOBAL, "a", encoding="utf-8") as f:
                f.write(message + "\n")
        except Exception as e:
            # Avoid recursive logging failure loops.
            print(f"[LOG_LISTENER_ERROR] Failed writing log message: {e}", file=sys.stderr)
            sys.stderr.flush()

def start_logging_listener(log_file_path):
    """Start queue-based centralized log writer for this process tree."""
    global LOG_QUEUE, LOG_LISTENER_THREAD, LOG_FILE_PATH_GLOBAL
    if LOG_QUEUE is not None and LOG_LISTENER_THREAD is not None and LOG_LISTENER_THREAD.is_alive():
        return True
    try:
        LOG_FILE_PATH_GLOBAL = log_file_path
        LOG_QUEUE = multiprocessing.Queue()
        LOG_LISTENER_THREAD = threading.Thread(target=_log_listener_worker, daemon=True)
        LOG_LISTENER_THREAD.start()
        return True
    except Exception as e:
        print(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Warning: Failed to start log listener, fallback to direct file writes. Error: {e}")
        LOG_QUEUE = None
        LOG_LISTENER_THREAD = None
        LOG_FILE_PATH_GLOBAL = None
        return False

def stop_logging_listener():
    """Gracefully stop centralized log writer thread."""
    global LOG_QUEUE, LOG_LISTENER_THREAD, LOG_FILE_PATH_GLOBAL
    try:
        if LOG_QUEUE is not None:
            LOG_QUEUE.put(None)
    except Exception:
        pass
    try:
        if LOG_LISTENER_THREAD is not None:
            LOG_LISTENER_THREAD.join(timeout=3)
    except Exception:
        pass
    LOG_QUEUE = None
    LOG_LISTENER_THREAD = None
    LOG_FILE_PATH_GLOBAL = None

def init_worker_logging(log_queue):
    """ProcessPool initializer: attach child process to main log queue."""
    global LOG_QUEUE
    LOG_QUEUE = log_queue

def setup_logging(log_file_path):
    log_dir = os.path.dirname(log_file_path)
    try:
        if log_dir: os.makedirs(log_dir, exist_ok=True)
        with open(log_file_path, "w", encoding="utf-8") as f:
            f.write(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Log file initialized.\n")
        return True
    except OSError as e:
        print(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] CRITICAL Error creating log directory or initializing log file {log_file_path}: {e}. Exiting.")
        return False

def log(message, log_file_path, console_output=True, level="INFO"):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lvl = str(level).upper() if level else "INFO"
    log_message = f"[{timestamp}] [{lvl}] {message}"
    if console_output:
        print(log_message, file=sys.stdout)
        sys.stdout.flush()

    # Preferred path: enqueue to a single writer (safe across processes).
    if LOG_QUEUE is not None:
        try:
            LOG_QUEUE.put(log_message)
            return
        except Exception:
            # Fallback to direct write if queue is unavailable.
            pass

    try:
        log_dir = os.path.dirname(log_file_path)
        if log_dir: os.makedirs(log_dir, exist_ok=True)
        with open(log_file_path, "a", encoding="utf-8") as f: f.write(log_message + "\n")
    except Exception as e:
        err_msg = f"[{timestamp}] CRITICAL: Failed to write to log file {log_file_path}: {e}"
        print(err_msg, file=sys.stderr)
        sys.stderr.flush()

def run_command(command, log_file_path, check=True):
    log(f"Running command: {command}", log_file_path, console_output=False, level="DEBUG")
    try:
        result = subprocess.run(['/bin/bash', '-c', command], check=check, capture_output=True, text=True, errors='ignore')
        if result.stdout: log(f"Stdout: {result.stdout.strip()}", log_file_path, console_output=False)
        if result.stderr: log(f"Stderr: {result.stderr.strip()}", log_file_path, console_output=False)
        return result
    except subprocess.CalledProcessError as e:
        log(f"Command failed (exit code {e.returncode}): {command}", log_file_path, level="ERROR")
        if e.stderr: log(f"Stderr: {e.stderr.strip()}", log_file_path, level="ERROR")
        return None
    except FileNotFoundError as e:
        log(f"Command or executable not found: {e}. Command: {command}", log_file_path, level="ERROR")
        return None
    except Exception as e:
        log(f"Unexpected exception while running command: {command}\nException: {e}", log_file_path, level="ERROR")
        return None

def count_sequences_fasta(fasta_path, log_file_path):
    if not os.path.exists(fasta_path) or os.path.getsize(fasta_path) == 0:
        return 0
    count_cmd = f"seqkit fx2tab --name --quiet {shlex.quote(fasta_path)} | wc -l"
    try:
        result = subprocess.run(['/bin/bash', '-c', count_cmd], check=True, capture_output=True, text=True)
        return int(result.stdout.strip())
    except (subprocess.CalledProcessError, ValueError, FileNotFoundError):
        qpath = shlex.quote(fasta_path)
        grep_cmd = f"zgrep -c '^>' {qpath}" if fasta_path.endswith(".gz") else f"grep -c '^>' {qpath}"
        try:
            result = subprocess.run(['/bin/bash', '-c', grep_cmd], check=True, capture_output=True, text=True)
            return int(result.stdout.strip())
        except (subprocess.CalledProcessError, ValueError, FileNotFoundError):
            log(f"Warning: Fallback grep count also failed for {fasta_path}. Returning 0.", log_file_path, console_output=False)
            return 0

def count_base_frequencies(fasta_file, log_file_path):
    sequence_count = 0
    valid_bases = ('A', 'T', 'C', 'G')
    open_func = gzip.open if fasta_file.endswith(".gz") else open
    read_mode = 'rt'
    analysis_limit = None
    try:
        if ':' in SEQKIT_RANGE:
            limit_str = SEQKIT_RANGE.split(':')[-1]
            if limit_str:
                analysis_limit = max(0, int(limit_str))
    except Exception:
        analysis_limit = None

    position_counters = []
    informative_totals = []
    ambiguous_char_counts = []
    n_char_counts = []
    unrecognized_char_counts = []

    def ensure_capacity(target_len):
        while len(position_counters) < target_len:
            position_counters.append(collections.Counter())
            informative_totals.append(0.0)
            ambiguous_char_counts.append(0)
            n_char_counts.append(0)
            unrecognized_char_counts.append(0)

    def process_sequence(sequence):
        if not sequence:
            return
        seq = sequence.upper()
        if analysis_limit is not None:
            seq = seq[:analysis_limit]
        if not seq:
            return
        ensure_capacity(len(seq))
        for idx, base in enumerate(seq):
            if base == 'N':
                n_char_counts[idx] += 1
                continue
            weights = IUPAC_BASE_WEIGHTS.get(base)
            if weights is None:
                unrecognized_char_counts[idx] += 1
                continue
            if base not in valid_bases:
                ambiguous_char_counts[idx] += 1
            informative_totals[idx] += 1.0
            for actual_base, weight in weights.items():
                position_counters[idx][actual_base] += weight

    try:
        with open_func(fasta_file, read_mode, errors='ignore') as f:
            current_seq_parts = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    sequence_count += 1
                    if current_seq_parts:
                        process_sequence("".join(current_seq_parts))
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line)
            if current_seq_parts:
                process_sequence("".join(current_seq_parts))
    except FileNotFoundError:
        log(f"Error: Fasta file not found at {fasta_file} for base frequency counting.", log_file_path)
        return None, 0
    except Exception as e:
        log(f"Error reading base frequencies from FASTA {fasta_file}: {e}", log_file_path)
        return None, 0

    analysis_length = len(position_counters)
    if analysis_length <= 0:
        return None, sequence_count

    freq_data = []
    for i in range(analysis_length):
        counter = position_counters[i]
        informative_total = informative_totals[i]
        if informative_total > 0:
            base_scores = {base: float(counter.get(base, 0.0)) for base in valid_bases}
            max_base = max(base_scores, key=base_scores.get)
            max_count = base_scores[max_base]
            max_freq_percent = (max_count / informative_total) * 100.0
        else:
            max_base = 'N'
            max_freq_percent = 0.0
        freq_data.append({
            'Pos': i + 1,
            'MaxBase': max_base,
            'MaxFreq': max_freq_percent,
            'InformativeReads': informative_total,
            'AmbiguousReads': ambiguous_char_counts[i],
            'NReads': n_char_counts[i],
            'UnrecognizedReads': unrecognized_char_counts[i]
        })

    if not freq_data:
        return None, sequence_count
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
        ax.set_xticks(range(0, int(plot_xmax) + 1, 20 if plot_xmax > 40 else 10))
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.set_xlabel("Position", fontsize=8)
        ax.set_ylabel("Max Base Frequency (%)", fontsize=8)
        ax.set_title(f"{entity_id} | {reads_in_entity:,} reads", fontsize=9, pad=10, fontweight='bold')

        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(
                handles,
                labels,
                loc='lower right',
                fontsize=6,
                frameon=False,
                handletextpad=0.2,
                labelspacing=0.3,
                borderpad=0.2,
            )

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
            os.path.join(input_dir, f)
            for f in os.listdir(input_dir)
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
    quoted_input_files = [shlex.quote(p) for p in input_files]
    use_zcat = input_files[0].endswith(".gz") if input_files else False
    cat_cmd = "zcat -f" if use_zcat else "cat"
    q_pool_out = shlex.quote(filtered_pool_file)
    q_range = shlex.quote(str(seq_range))
    q_maxee = shlex.quote(str(max_ee))

    filter_cmd = (
        f"{cat_cmd} {' '.join(quoted_input_files)} | "
        f"seqkit subseq --quiet -r {q_range} | "
        f"vsearch --fastq_filter - --fastq_maxee {q_maxee} --fastaout {q_pool_out} --quiet"
    )

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

    cluster_cmd = (
        f"vsearch --cluster_fast {shlex.quote(filtered_pool_file)} --id {cluster_id_threshold} "
        f"--centroids {shlex.quote(centroids_file)} --clusters {shlex.quote(cluster_file_prefix)} "
        f"--threads {num_threads} --quiet"
    )

    cluster_result = run_command(cluster_cmd, log_file_path)
    if cluster_result is None or cluster_result.returncode != 0:
        log(f"Error: VSEARCH clustering failed. Exiting.", log_file_path); return None

    base_prefix = os.path.basename(cluster_file_prefix)
    glob_pattern = os.path.join(clusters_dir, f"{base_prefix}*")
    initial_cluster_files = [
        f for f in glob.glob(glob_pattern)
        if os.path.basename(f).replace(base_prefix, "").replace(".fasta", "").replace(".fa", "").isdigit()
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
            try:
                os.remove(cluster_file_path)
            except OSError:
                pass
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
    q_in = shlex.quote(sample_fastq_path)
    q_out = shlex.quote(output_fasta)
    q_range = shlex.quote(str(seq_range))
    q_maxee = shlex.quote(str(max_ee))
    filter_cmd = (
        f"{cat_cmd} {q_in} | "
        f"seqkit subseq --quiet -r {q_range} | "
        f"vsearch --fastq_filter - --fastq_maxee {q_maxee} --fastaout {q_out} --quiet"
    )
    filter_result = run_command(filter_cmd, log_file_path)

    if filter_result is None or filter_result.returncode != 0:
        log(f"Error: Filtering command failed for sample '{sample_id}'.", log_file_path)
        return None, 0

    read_count = count_sequences_fasta(output_fasta, log_file_path)
    if read_count == 0:
        try:
            os.remove(output_fasta)
        except OSError:
            pass
        return None, 0
    return output_fasta, read_count

def analyze_entities_parallel(iterable_of_args, process_function_partial, entity_type_desc, num_workers, log_file_path):
    num_entities_to_analyze = len(iterable_of_args)
    log(f"Parallel Analysis: Submitting {num_entities_to_analyze} {entity_type_desc.lower()} tasks ({num_workers} workers)...", log_file_path)

    analysis_results_dict = {}
    futures_list = []

    executor_kwargs = {'max_workers': num_workers}
    if LOG_QUEUE is not None:
        executor_kwargs['initializer'] = init_worker_logging
        executor_kwargs['initargs'] = (LOG_QUEUE,)

    with ProcessPoolExecutor(**executor_kwargs) as executor:
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

    if not analysis_results_list and num_entities_to_analyze > 0:
        log(f"Error: No {entity_type_desc.lower()} successfully analyzed. Exiting.", log_file_path); return None, []

    entity_consensus_results = [(res['id'], res['consensus_seq']) for res in analysis_results_list if res.get('consensus_seq')]
    log(f"Parallel Analysis: Generated {len(entity_consensus_results)} non-empty consensus sequences.", log_file_path)
    return analysis_results_list, entity_consensus_results

def generate_plots_for_entities(analysis_results_list, plots_dir, log_file_path, entity_type_desc):
    if not analysis_results_list:
        return 0, 0

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
            else:
                plots_failed_count += 1
        else:
            plots_failed_count += 1

        if (i + 1) % LOG_INCREMENT_LOOPS == 0 or (i + 1) == num_results:
            log(f"Plot Generation: Processed {i+1}/{num_results} plots...", log_file_path, console_output=False)

    log(f"Plot Generation: Finished. Generated {plots_generated_count}, Failed/Skipped {plots_failed_count} for {entity_type_desc.lower()}.", log_file_path)
    return plots_generated_count, plots_failed_count

def run_blast_pipeline(entity_consensus_results, blast_input_fasta_name, blast_db_path, blast_evalue_param, blast_outfmt_param, blast_output_tsv_path, best_hits_csv_path, num_threads_param, log_file_path, entity_type_desc, analysis_results_list=None, blast_db_meta=None, sample_suffix=None):
    if not entity_consensus_results:
        log(f"BLAST: Skipping for {entity_type_desc} - No consensus sequences.", log_file_path)
        return False, None

    blast_dir = os.path.dirname(blast_output_tsv_path)
    os.makedirs(blast_dir, exist_ok=True)
    blast_input_fasta_path = os.path.join(blast_dir, blast_input_fasta_name)

    try:
        with open(blast_input_fasta_path, "w") as f_out:
            for entity_id, sequence in entity_consensus_results:
                f_out.write(f">{entity_id}\n{sequence}\n")
    except IOError as e:
        log(f"Error (BLAST): Error writing BLAST input {blast_input_fasta_path}: {e}. Skipping BLAST.", log_file_path)
        return False, None

    meta = blast_db_meta or {}
    db_label = meta.get('label', os.path.basename(blast_db_path))
    read_orientation = infer_read_orientation(sample_suffix)
    max_len = max(len(s) for _, s in entity_consensus_results) if entity_consensus_results else 0
    task = "blastn-short" if 0 < max_len < 50 else "blastn"
    log(
        f"BLAST: {len(entity_consensus_results)} {entity_type_desc.lower()} | db={db_label!r} | {task} | "
        f"e<{blast_evalue_param} | word_size={BLAST_WORD_SIZE} | strand=both | "
        f"read={read_orientation or 'unknown'} | orientation_filter={'on' if read_orientation else 'off'}",
        log_file_path,
    )
    log(f"BLAST: database path {blast_db_path}", log_file_path, console_output=False)
    cmd = (
        f"blastn -query {shlex.quote(blast_input_fasta_path)} -db {shlex.quote(blast_db_path)} "
        f"-out {shlex.quote(blast_output_tsv_path)} -outfmt {shlex.quote(blast_outfmt_param)} "
        f"-evalue {shlex.quote(str(blast_evalue_param))} -task {shlex.quote(task)} "
        f"-word_size {shlex.quote(str(BLAST_WORD_SIZE))} -num_threads {shlex.quote(str(num_threads_param))} -strand both"
    )

    result = run_command(cmd, log_file_path)
    if result is None or result.returncode != 0:
        log(f"Warning (BLAST): BLAST command may have failed. Output will be processed if available.", log_file_path)

    log(f"BLAST: parse TSV -> CSV ({entity_type_desc})", log_file_path)
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
        return False, None

    primer_call_lines = None

    all_queries_df = pd.DataFrame(all_submitted_queries_data)

    hits_from_blast_df = pd.DataFrame()
    ranked_hits_df = pd.DataFrame()
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
                        # Prioritize hits per query (primary: bit score; tie-breakers: identity, length, evalue).
                        ranked_hits_df = hits_from_blast_df.copy()
                        numeric_cols = ['bit_score', '%_identity', 'alignment_length', 'evalue']
                        for col in numeric_cols:
                            if col in ranked_hits_df.columns:
                                ranked_hits_df[col] = pd.to_numeric(ranked_hits_df[col], errors='coerce')
                        sort_cols = [c for c in ['query_id', 'bit_score', '%_identity', 'alignment_length', 'evalue'] if c in ranked_hits_df.columns]
                        sort_ascending = []
                        for col in sort_cols:
                            if col == 'query_id':
                                sort_ascending.append(True)
                            elif col == 'evalue':
                                sort_ascending.append(True)
                            else:
                                sort_ascending.append(False)
                        ranked_hits_df = ranked_hits_df.sort_values(sort_cols, ascending=sort_ascending, na_position='last')
                        ranked_hits_df['hit_rank'] = ranked_hits_df.groupby('query_id', sort=False).cumcount() + 1
                        ranked_hits_df['is_top_suspect'] = ranked_hits_df['hit_rank'] == 1
                        ranked_hits_df = filter_ranked_hits_by_orientation(
                            ranked_hits_df, read_orientation, db_label, log_file_path
                        )
                else:
                    missing_cols = [col for col in ['bit_score', 'query_id'] if col not in df_raw_blast.columns]
                    log(f"Warning (BLAST): Key column(s) {missing_cols} not found in raw BLAST output. Cannot process hits.", log_file_path, level="WARNING")
        except Exception as e:
            log(f"Error (BLAST): Failed to process raw BLAST output {os.path.basename(blast_output_tsv_path)}: {e}", log_file_path, level="ERROR")

    if not all_queries_df.empty:
        if not ranked_hits_df.empty:
            final_df = pd.merge(all_queries_df, ranked_hits_df, on="query_id", how="left")
        else:
            final_df = all_queries_df.copy()
            for col_name in csv_column_names:
                if col_name not in final_df.columns:
                    final_df[col_name] = pd.NA
            final_df['hit_rank'] = pd.NA
            final_df['is_top_suspect'] = pd.NA

        pfx = meta.get('prefix', blast_db_path)
        lbl = meta.get('label', os.path.basename(str(pfx)))
        res = meta.get('resolution', 'N/A')
        final_df['blast_db_prefix'] = pfx
        final_df['blast_db_label'] = lbl
        final_df['blast_db_resolution'] = res
        final_df['top_suspect_recommendation'] = final_df['is_top_suspect'].map({True: 'Y', False: 'N'})
        final_df.loc[final_df['is_top_suspect'].isna(), 'top_suspect_recommendation'] = 'N/A'

        if '%_identity' in final_df.columns:
            final_df['%_identity'] = pd.to_numeric(final_df['%_identity'], errors='coerce').round(2)
        if 'bit_score' in final_df.columns:
            final_df['bit_score'] = pd.to_numeric(final_df['bit_score'], errors='coerce').round(2)

        # Export compact top-suspect table (one recommended hit per query).
        top_suspect_output_path = best_hits_csv_path.replace("_final_blast_results.csv", "_top_suspect.csv")
        try:
            top_cols = ['query_id', 'query_reads_count', 'subject_id', '%_identity', 'alignment_length', 'evalue', 'bit_score']
            available_top_cols = [c for c in top_cols if c in final_df.columns]
            top_rows = final_df[final_df['is_top_suspect'] == True][available_top_cols].copy()  # noqa: E712
            top_suspect_df = pd.merge(all_queries_df, top_rows, on=['query_id', 'query_reads_count'], how='left')
            top_suspect_df['blast_db_prefix'] = pfx
            top_suspect_df['blast_db_label'] = lbl
            top_suspect_df['blast_db_resolution'] = res
            top_suspect_df['status'] = top_suspect_df['subject_id'].isna().map({True: 'miss', False: 'hit'})
            top_suspect_df = top_suspect_df.rename(
                columns={k: v for k, v in BLAST_EXPORT_RENAME.items() if k in top_suspect_df.columns}
            )
            for col in TOP_SUSPECT_EXPORT_COLS:
                if col not in top_suspect_df.columns:
                    top_suspect_df[col] = pd.NA
            top_suspect_df = top_suspect_df[list(TOP_SUSPECT_EXPORT_COLS)].astype(object)
            top_suspect_df = top_suspect_df.where(pd.notna(top_suspect_df), '')
            top_suspect_df.to_csv(top_suspect_output_path, index=False)
            log(f"BLAST: summary -> {os.path.basename(top_suspect_output_path)} ({len(top_suspect_df)} rows)", log_file_path)
            primer_call_lines = write_primer_call_summary(
                os.path.dirname(best_hits_csv_path),
                top_suspect_df,
                entity_type_desc,
                lbl,
                read_orientation,
                log_file_path,
            )
        except Exception as e:
            log(f"Warning (BLAST): Failed to save top-suspect summary: {e}", log_file_path, level="WARNING")

        for col in final_df.columns:
            final_df[col] = final_df[col].astype("object").fillna("N/A").infer_objects()

        desired_column_order = ['blast_db_prefix', 'blast_db_label', 'blast_db_resolution', 'query_id', 'query_reads_count', 'hit_rank', 'top_suspect_recommendation']
        for mapped_col in csv_column_names:
            if mapped_col != 'query_id' and mapped_col not in desired_column_order:
                desired_column_order.append(mapped_col)

        for col in desired_column_order:
            if col not in final_df.columns:
                final_df[col] = "N/A"
        final_df = final_df[desired_column_order]

        try:
            export_df = prepare_blast_results_for_export(final_df)
            export_df.to_csv(best_hits_csv_path, index=False)
            n_queries = export_df['query'].nunique() if 'query' in export_df.columns else 0
            log(
                f"BLAST: hits table -> {os.path.basename(best_hits_csv_path)} ({len(export_df)} rows, {n_queries} queries)",
                log_file_path,
            )
            processed_successfully = True
        except Exception as e:
            log(f"Error (BLAST): Failed to save final BLAST report to {os.path.basename(best_hits_csv_path)}: {e}", log_file_path)
    else:
        log(f"BLAST: No submitted queries to generate a final report for {entity_type_desc}.", log_file_path)

    return processed_successfully, primer_call_lines

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
    blast_db_group = parser.add_mutually_exclusive_group()
    blast_db_group.add_argument(
        "--blast-db",
        dest="blast_db",
        default=None,
        choices=["its", "16s"],
        type=_blast_cli_builtin_db,
        metavar="{its,16s}",
        help="Built-in primer library: its or 16s. Mutually exclusive with --blast-db-dir.",
    )
    blast_db_group.add_argument(
        "--blast-db-dir",
        dest="blast_db_dir",
        default=None,
        metavar="DIR",
        help="Directory containing a custom BLAST DB (auto-detects a single prefix). Mutually exclusive with --blast-db.",
    )
    args = parser.parse_args()

    input_dir, output_dir, sample_suffix = args.input_dir, args.output_dir, args.sample_suffix
    THREADS = args.threads
    BLAST_EVALUE = args.e_value
    run_pooled_mode = args.pool

    pipeline_start_time = time.time()
    logs_dir = os.path.join(output_dir, "logs")
    log_file = os.path.join(logs_dir, "pipeline.log")
    if not setup_logging(log_file): sys.exit(1)
    start_logging_listener(log_file)
    atexit.register(stop_logging_listener)

    blast_db_to_use, blast_db_meta = resolve_blast_database(args)
    if blast_db_meta.get('error'):
        log(f"CRITICAL: {blast_db_meta['error']}", log_file)
        sys.exit(1)

    log("=== APHunter start ===", log_file)
    log(
        f"In {os.path.abspath(input_dir)} -> Out {os.path.abspath(output_dir)} | "
        f"{'pool' if run_pooled_mode else 'per-sample'} | threads={THREADS}",
        log_file,
    )
    db_label_for_filename = sanitize_db_label_for_filename(blast_db_meta.get('label', ''))
    log(
        f"BLAST db={blast_db_meta.get('label', '')!r} from={blast_db_meta.get('resolution', '')} | {blast_db_to_use}",
        log_file,
    )
    log(f"BLAST settings: e-value <= {BLAST_EVALUE} | word_size={BLAST_WORD_SIZE} | strand=both", log_file, console_output=False)
    read_orientation = infer_read_orientation(sample_suffix)
    if read_orientation:
        log(f"BLAST orientation filter: {read_orientation} read -> keep {read_orientation} primers only", log_file, console_output=False)
    else:
        log(f"BLAST orientation filter: skipped (could not infer read orientation from suffix {sample_suffix!r})", log_file, level="WARNING")

    final_dir = os.path.join(output_dir, OUT_FINAL)
    if run_pooled_mode:
        pool_work_dir = os.path.join(output_dir, OUT_POOL_WORK)
        clusters_dir = os.path.join(output_dir, OUT_POOL_CLUSTERS)
        stats_dir = os.path.join(output_dir, OUT_POOL_STATS)
        plots_dir = os.path.join(output_dir, OUT_POOL_PLOTS)
        blast_work_dir = os.path.join(output_dir, OUT_POOL_BLAST_WORK)
        all_subdirs = [pool_work_dir, clusters_dir, stats_dir, plots_dir, blast_work_dir, final_dir, logs_dir]
        filtered_samples_output_dir = None
    else:
        filtered_samples_output_dir = os.path.join(output_dir, OUT_NONPOOL_FILTERED)
        stats_dir = os.path.join(output_dir, OUT_NONPOOL_STATS)
        plots_dir = os.path.join(output_dir, OUT_NONPOOL_PLOTS)
        blast_work_dir = os.path.join(output_dir, OUT_NONPOOL_BLAST_WORK)
        clusters_dir = None
        all_subdirs = [filtered_samples_output_dir, stats_dir, plots_dir, blast_work_dir, final_dir, logs_dir]

    log(
        f"Output layout: numbered folders = intermediates; {OUT_FINAL}/ = BLAST result tables.",
        log_file,
    )

    analysis_results_list = []
    entity_consensus_results = []
    final_blast_output_csv = ""
    primer_call_lines = None
    initial_sample_count = 0
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
        filtered_pool_file = os.path.join(pool_work_dir, "pooled_filtered.fasta")
        all_centroids_file = os.path.join(clusters_dir, "all_centroids.fasta")
        cluster_file_prefix = os.path.join(clusters_dir, "cluster.")

        log(f"Step {step_num}: Pooling & Filtering reads...", log_file)
        filtered_pool_file, pooled_read_count = run_pooling_and_filtering_step(input_files, filtered_pool_file, SEQKIT_RANGE, MAX_EE, log_file)
        if filtered_pool_file is None: sys.exit(1)
        step_num += 1

        min_reads_cluster = max(MIN_READS_FLOOR, pooled_read_count // 2000 if pooled_read_count > 0 else 1)

        log(f"Step {step_num}: Clustering reads (ID >= {CLUSTER_IDENTITY_THRESHOLD})...", log_file)
        initial_cluster_files = run_clustering_step(filtered_pool_file, CLUSTER_IDENTITY_THRESHOLD, all_centroids_file, cluster_file_prefix, THREADS, log_file)
        if initial_cluster_files is None: sys.exit(1)
        try:
            os.remove(filtered_pool_file)
        except OSError:
            pass
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
            blast_output_tsv_pooled = os.path.join(blast_work_dir, f"pooled_{db_label_for_filename}_blast_output.tsv")
            final_blast_output_csv = os.path.join(final_dir, f"pooled_{db_label_for_filename}_final_blast_results.csv")
            _, primer_call_lines = run_blast_pipeline(
                entity_consensus_results,
                f"cluster_consensus_for_blast_{db_label_for_filename}.fasta",
                blast_db_to_use,
                BLAST_EVALUE,
                BLAST_OUTFMT,
                blast_output_tsv_pooled,
                final_blast_output_csv,
                THREADS,
                log_file,
                "Clusters",
                analysis_results_list,
                blast_db_meta=blast_db_meta,
                sample_suffix=sample_suffix,
            )

    else:  # NON-POOLED (PER-SAMPLE) MODE
        log(f"--- Starting NON-POOLED (Per-Sample) MODE Pipeline ---", log_file)

        log(f"Step {step_num}: Filtering individual samples (Min reads: {MIN_READS_PER_SAMPLE_NON_POOL})...", log_file)
        filtered_sample_fastas_args = []
        samples_discarded_count = 0
        for i, fq_path in enumerate(input_files):
            sample_id = get_sample_id_from_path(fq_path, sample_suffix)
            filtered_fasta, read_count = filter_single_sample_reads(fq_path, sample_id, filtered_samples_output_dir, SEQKIT_RANGE, MAX_EE, log_file)
            if filtered_fasta and read_count >= MIN_READS_PER_SAMPLE_NON_POOL:
                filtered_sample_fastas_args.append((filtered_fasta, sample_id))
            else:
                samples_discarded_count += 1

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
        blast_output_tsv_nonpool = os.path.join(blast_work_dir, f"per_sample_{db_label_for_filename}_blast_output.tsv")
        final_blast_output_csv = os.path.join(final_dir, f"per_sample_{db_label_for_filename}_final_blast_results.csv")
        _, primer_call_lines = run_blast_pipeline(
            entity_consensus_results,
            f"sample_consensus_for_blast_{db_label_for_filename}.fasta",
            blast_db_to_use,
            BLAST_EVALUE,
            BLAST_OUTFMT,
            blast_output_tsv_nonpool,
            final_blast_output_csv,
            THREADS,
            log_file,
            "Samples",
            analysis_results_list,
            blast_db_meta=blast_db_meta,
            sample_suffix=sample_suffix,
        )

    elapsed_seconds = time.time() - pipeline_start_time
    elapsed_minutes, elapsed_sec_rem = divmod(elapsed_seconds, 60)
    log("=== Summary ===", log_file)
    log(
        f"Time {int(elapsed_minutes)}m {elapsed_sec_rem:.1f}s | out={os.path.abspath(output_dir)} | "
        f"inputs={initial_sample_count} | kept={final_entities_count} | analyzed={final_entities_analyzed_count} | "
        f"consensus_BLAST={final_entities_for_blast_count}",
        log_file,
    )
    log(f"Primary BLAST CSVs: {os.path.abspath(final_dir)}/", log_file)

    if entity_consensus_results and os.path.exists(final_blast_output_csv) and os.path.getsize(final_blast_output_csv) > 0:
        try:
            with open(final_blast_output_csv, 'r', newline='') as f_blast_csv:
                header_line = f_blast_csv.readline()
                actual_headers_in_file = next(csv.reader([header_line]))
                reads_ok = ('reads' in actual_headers_in_file) or ('query_reads_count' in actual_headers_in_file)
                if reads_ok:
                    log(f"BLAST CSV header OK ({os.path.basename(final_blast_output_csv)}).", log_file, console_output=False)
                else:
                    log(f"Warning: expected reads column missing in {os.path.basename(final_blast_output_csv)}.", log_file)
                meta_ok = (
                    all(h in actual_headers_in_file for h in ('db_path', 'db', 'db_from'))
                    or all(h in actual_headers_in_file for h in ('blast_db_prefix', 'blast_db_label', 'blast_db_resolution'))
                )
                if meta_ok:
                    log(f"BLAST DB metadata columns present.", log_file, console_output=False)

            with open(final_blast_output_csv, 'r') as f_blast_csv_count:
                num_total_entries = max(0, sum(1 for _ in f_blast_csv_count) - 1)
            log(f"BLAST table rows: {num_total_entries} ({os.path.basename(final_blast_output_csv)})", log_file)
        except Exception as e_count:
            log(f"Final BLAST Report: {os.path.basename(final_blast_output_csv)} (exists, but count/header check failed: {e_count})", log_file)
    elif entity_consensus_results:
        log(f"BLAST was run, but final report '{os.path.basename(final_blast_output_csv)}' is missing or empty.", log_file)
    else:
        log(f"BLAST step was skipped (no consensus sequences).", log_file)

    log(f"log file: {log_file}", log_file)
    if primer_call_lines:
        log("=== Primer call ===", log_file)
        for line in primer_call_lines:
            log(line, log_file)
    log("=== end ===", log_file)

if __name__ == "__main__":
    main()
