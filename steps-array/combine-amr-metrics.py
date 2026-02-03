#!/usr/bin/env python3

import argparse
import os
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd


BOWTIE_TOTAL_READS_RE = re.compile(r"^(\d+)\s+reads; of these:", re.MULTILINE)
BOWTIE_OVERALL_ALN_RE = re.compile(r"([0-9.]+)%\s+overall alignment rate", re.MULTILINE)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def find_files(root: Path, pattern: str) -> List[Path]:
    return sorted(root.rglob(pattern))


def sample_from_amr(path: Path) -> str:
    # amrfp_<SAMPLE>.tsv
    return re.sub(r"^amrfp_", "", path.stem)


def sample_from_cov(path: Path) -> str:
    # covStats_<SAMPLE>.txt
    return re.sub(r"^covStats_", "", path.stem)


def sample_from_maplog(path: Path) -> str:
    # bowtie2_mapLog_<SAMPLE>.log
    return re.sub(r"^bowtie2_mapLog_", "", path.stem)


def sample_from_spf_pe(path: Path) -> str:
    # <SAMPLE>_singlem.spf.tsv
    return re.sub(r"_singlem\.pe_spf\.tsv$", "", path.name)

def sample_from_spf_se(path: Path) -> str:
    # <SAMPLE>_singlem.spf.tsv
    return re.sub(r"_singlem\.se_spf\.tsv$", "", path.name)


def read_amr_file(path: Path, amr_cols: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Read AMRFinderPlus TSV. If empty/unreadable, return a 1-row NA frame with known columns.
    Adds a 'Name' column with sample ID.
    """
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        df = pd.DataFrame()
    
    if df.empty:
        if amr_cols is None:
            # If we can't infer columns, return minimal stub
            return pd.DataFrame({"Name": [sample_from_amr(path)]})
        stub = {c: [pd.NA] for c in amr_cols}
        df = pd.DataFrame(stub)
    
    df["Name"] = sample_from_amr(path)
    return df


def parse_bowtie2_log(path: Path) -> Tuple[Optional[int], Optional[float]]:
    """Extract total reads and overall alignment percent from bowtie2 log."""
    try:
        txt = path.read_text()
    except Exception:
        return None, None
    
    m1 = BOWTIE_TOTAL_READS_RE.search(txt)
    m2 = BOWTIE_OVERALL_ALN_RE.search(txt)
    
    total_reads = int(m1.group(1)) if m1 else None
    overall_aln_percent = float(m2.group(1)) if m2 else None
    return total_reads, overall_aln_percent


def read_covstats_file(path: Path, sample: str) -> pd.DataFrame:
    """
    Read samtools coverage output. Handles '#rname' header by stripping '#'.
    Normalises first column to 'rname'. Keeps the key numeric fields if present.
    """
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return pd.DataFrame()
    
    if df.empty:
        return df
    
    # Strip leading '#'
    df.columns = [c.lstrip("#") for c in df.columns]
    
    # Ensure rname exists
    if "rname" not in df.columns:
        # often the first column is rname
        df = df.rename(columns={df.columns[0]: "rname"})
    
    keep = ["rname", "numreads", "covbases", "coverage", "meandepth"]
    present = [c for c in keep if c in df.columns]
    df = df[present].copy()
    df["sample"] = sample
    return df


def read_spf_file_pe(path: Path, sample: str) -> pd.DataFrame:
    """
    Read singlem prokaryotic_fraction output TSV and extract bacterial_archaeal_bases and metagenome_size.
    Returns both as separate columns with '_pe' suffix.
    """
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return pd.DataFrame({
            "sample": [sample],
            "bacterial_archaeal_bases_pe": [pd.NA],
            "metagenome_size_pe": [pd.NA]
        })
    
    if df.empty:
        return pd.DataFrame({
            "sample": [sample],
            "bacterial_archaeal_bases_pe": [pd.NA],
            "metagenome_size_pe": [pd.NA]
        })
    
    bacterial_col = "bacterial_archaeal_bases" if "bacterial_archaeal_bases" in df.columns else None
    metagenome_col = "metagenome_size" if "metagenome_size" in df.columns else None
    
    bacterial_val = df.iloc[0][bacterial_col] if bacterial_col else pd.NA
    metagenome_val = df.iloc[0][metagenome_col] if metagenome_col else pd.NA
    
    return pd.DataFrame({
        "sample": [sample],
        "bacterial_archaeal_bases_pe": [bacterial_val],
        "metagenome_size_pe": [metagenome_val]
    })


def read_spf_file_s1(path: Path, sample: str) -> pd.DataFrame:
    """
    Read singlem prokaryotic_fraction output TSV and extract bacterial_archaeal_bases and metagenome_size.
    Returns both as separate columns with '_s1' suffix.
    """
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return pd.DataFrame({
            "sample": [sample],
            "bacterial_archaeal_bases_s1": [pd.NA],
            "metagenome_size_s1": [pd.NA]
        })
    
    if df.empty:
        return pd.DataFrame({
            "sample": [sample],
            "bacterial_archaeal_bases_s1": [pd.NA],
            "metagenome_size_s1": [pd.NA]
        })
    
    bacterial_col = "bacterial_archaeal_bases" if "bacterial_archaeal_bases" in df.columns else None
    metagenome_col = "metagenome_size" if "metagenome_size" in df.columns else None
    
    bacterial_val = df.iloc[0][bacterial_col] if bacterial_col else pd.NA
    metagenome_val = df.iloc[0][metagenome_col] if metagenome_col else pd.NA
    
    return pd.DataFrame({
        "sample": [sample],
        "bacterial_archaeal_bases_s1": [bacterial_val],
        "metagenome_size_s1": [metagenome_val]
    })

def read_spf_file_s2(path: Path, sample: str) -> pd.DataFrame:
    """
    Read singlem prokaryotic_fraction output TSV and extract bacterial_archaeal_bases and metagenome_size.
    Returns both as separate columns with '_s2' suffix.
    """
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return pd.DataFrame({
            "sample": [sample],
            "bacterial_archaeal_bases_s2": [pd.NA],
            "metagenome_size_s2": [pd.NA]
        })
    
    if df.empty:
        return pd.DataFrame({
            "sample": [sample],
            "bacterial_archaeal_bases_s2": [pd.NA],
            "metagenome_size_s2": [pd.NA]
        })
    
    bacterial_col = "bacterial_archaeal_bases" if "bacterial_archaeal_bases" in df.columns else None
    metagenome_col = "metagenome_size" if "metagenome_size" in df.columns else None
    
    bacterial_val = df.iloc[0][bacterial_col] if bacterial_col else pd.NA
    metagenome_val = df.iloc[0][metagenome_col] if metagenome_col else pd.NA
    
    return pd.DataFrame({
        "sample": [sample],
        "bacterial_archaeal_bases_s2": [bacterial_val],
        "metagenome_size_s2": [metagenome_val]
    })
    

def to_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Combine AMRFinderPlus + bowtie2 + coverage (+ optional SingleM SPF) into normalised AMR metrics."
    )
    parser.add_argument(
        "--results-root",
        default=None,
        help="Root directory to search recursively (default: $OUT_DIR or <submit_dir>/results).",
    )
    parser.add_argument(
        "--out-csv",
        default=None,
        help="Output CSV path (default: <results-root>/amr-reads-summary.csv).",
    )
    parser.add_argument(
        "--singlem",
        choices=["auto", "on", "off"],
        default="auto",
        help=(
            "Control use of SingleM SPF output. "
            "'auto' uses SPF if files exist, 'on' requires SPF files, 'off' ignores SPF."
        ),
    )
    args = parser.parse_args()
    
    submit_dir = Path(os.environ.get("SLURM_SUBMIT_DIR", os.getcwd())).resolve()
    if args.results_root:
        root = Path(args.results_root).resolve()
    else:
        out_dir = os.environ.get("OUT_DIR", str(submit_dir / "results"))
        root = Path(out_dir).resolve()
    
    if args.out_csv:
        out_csv = Path(args.out_csv).resolve()
    else:
        out_csv = root / "amr-reads-summary.csv"
    
    if not root.exists():
        eprint(f"ERROR: results-root does not exist: {root}")
        return 2
    
    # Discover files
    amr_files = find_files(root, "amrfp_*.tsv")
    cov_files = find_files(root, "covStats_*.txt")
    map_files = find_files(root, "bowtie2_mapLog_*.log")
    spf_files_pe = find_files(root, "*_singlem.pe_spf.tsv")
    spf_files_s1 = find_files(root, "*_singlem.s1_spf.tsv")
    spf_files_s2 = find_files(root, "*_singlem.s2_spf.tsv")
    
    # Apply SingleM control
    if args.singlem == "auto":
        use_singlem = bool(spf_files_pe or spf_files_s2 or spf_files_s2)
    else:
        use_singlem = (args.singlem == "on")
    
    require_singlem = (args.singlem == "on")
    
    if require_singlem:
        if not spf_files_pe:
            eprint("ERROR: --singlem on was set but no *_singlem.pe_spf.tsv files were found.")
            return 2
        if not spf_files_s1:
            eprint("ERROR: --singlem on was set but no *_singlem.s1_spf.tsv files were found.")
            return 2
        if not spf_files_s2:
            eprint("ERROR: --singlem on was set but no *_singlem.s2_spf.tsv files were found.")
            return 2
    elif not use_singlem:
        spf_files_pe = []
        spf_files_se = []
    
    
    eprint(f"Found AMR files:     {len(amr_files)}")
    eprint(f"Found covStats:      {len(cov_files)}")
    eprint(f"Found map logs:      {len(map_files)}")
    eprint(f"Found SPF files:     {len(spf_files_pe)} (mode={args.singlem})")
    eprint(f"Found SPF files:     {len(spf_files_s1)} (mode={args.singlem})")
    eprint(f"Found SPF files:     {len(spf_files_s2)} (mode={args.singlem})")
    
    if not amr_files:
        eprint("ERROR: No amrFP_out_*.tsv files found.")
        return 2
    
    # Infer AMR columns from first readable file
    amr_cols = None
    for f in amr_files:
        try:
            tmp = pd.read_csv(f, sep="\t", nrows=1, dtype=str)
            if not tmp.empty:
                amr_cols = list(tmp.columns)
                break
        except Exception:
            continue
    
    # Read AMR
    amr_list = [read_amr_file(f, amr_cols=amr_cols) for f in amr_files]
    df_amr = pd.concat(amr_list, ignore_index=True)
    
    df_amr.rename(columns={
        "Contig id": "Contig.id",
        "Alignment length": "Alignment.length",
        "Reference sequence length": "Reference.sequence.length"
    }, inplace=True)
    
    if "Contig.id" not in df_amr.columns:
        eprint("ERROR: AMR table does not contain 'Contig.id'. Check AMRFinderPlus output format.")
        return 2
        
    # Bowtie2 logs -> per-sample totals
    bow_rows = []
    for f in map_files:
        sample = sample_from_maplog(f)
        total_reads, overall_aln_percent = parse_bowtie2_log(f)
        bow_rows.append(
            {"sample": sample, "total_reads": total_reads, "overall_aln_percent": overall_aln_percent}
        )
    df_bow = pd.DataFrame(bow_rows)
    
    # Coverage -> per contig
    cov_rows = []
    for f in cov_files:
        sample = sample_from_cov(f)
        cov_df = read_covstats_file(f, sample)
        if not cov_df.empty:
            cov_rows.append(cov_df)
    df_cov = pd.concat(cov_rows, ignore_index=True) if cov_rows else pd.DataFrame()
    
    # Merge AMR + bowtie
    df_amr = df_amr.merge(df_bow, how="left", left_on="Name", right_on="sample").drop(
        columns=["sample"], errors="ignore"
    )
    
    # Merge SPF_pe if enabled and present
    if spf_files_pe:
        spf_rows = []
        for f in spf_files_pe:
            sample = sample_from_spf_pe(f)
            spf_rows.append(read_spf_file_pe(f, sample))
        df_spf_pe = pd.concat(spf_rows, ignore_index=True) if spf_rows else pd.DataFrame()
    
        df_amr = df_amr.merge(df_spf_pe, how="left", left_on="Name", right_on="sample").drop(
            columns=["sample"], errors="ignore"
        )
    else:
        df_amr["bacterial_archaeal_bases_pe"] = pd.NA
        df_amr["metagenome_size_pe"] = pd.NA
    
    # Merge SPF_s1 if enabled and present
    if spf_files_s1:
        spf_rows = []
        for f in spf_files_s1:
            sample = sample_from_spf_s1(f)
            spf_rows.append(read_spf_file_s1(f, sample))
        df_spf_s1 = pd.concat(spf_rows, ignore_index=True) if spf_rows else pd.DataFrame()
    
        df_amr = df_amr.merge(df_spf_s1, how="left", left_on="Name", right_on="sample").drop(
            columns=["sample"], errors="ignore"
        )
    else:
        df_amr["bacterial_archaeal_bases_s1"] = pd.NA
        df_amr["metagenome_size_s1"] = pd.NA
    
    # Merge SPF_s2 if enabled and present
    if spf_files_s2:
        spf_rows = []
        for f in spf_files_s2:
            sample = sample_from_spf_s2(f)
            spf_rows.append(read_spf_file_s2(f, sample))
        df_spf_s2 = pd.concat(spf_rows, ignore_index=True) if spf_rows else pd.DataFrame()
    
        df_amr = df_amr.merge(df_spf_s2, how="left", left_on="Name", right_on="sample").drop(
            columns=["sample"], errors="ignore"
        )
    else:
        df_amr["bacterial_archaeal_bases_s2"] = pd.NA
        df_amr["metagenome_size_s2"] = pd.NA
        
    # Merge coverage on (Name, Contig.id) ↔ (sample, rname)
    if not df_cov.empty:
        df_amr = df_amr.merge(
            df_cov,
            how="left",
            left_on=["Name", "Contig.id"],
            right_on=["sample", "rname"],
        ).drop(columns=["sample", "rname"], errors="ignore")
    else:
        for c in ["numreads", "covbases", "coverage", "meandepth"]:
            if c not in df_amr.columns:
                df_amr[c] = pd.NA
    
    # Numeric conversions
    for c in [
        "numreads",
        "covbases",
        "coverage",
        "meandepth",
        "Alignment.length",
        "Reference.sequence.length",
        "total_reads",
        "overall_aln_percent",
        "prokaryotic_fraction",
    ]:
        if c in df_amr.columns:
            df_amr[c] = to_numeric(df_amr[c])
    
    # Metrics
    df_amr["base_cov_amr_per_million_reads"] = (
        df_amr["meandepth"] * df_amr["Alignment.length"] / (df_amr["total_reads"] / 1e6)
    )
    
    df_amr["rpkm"] = 1e9 * df_amr["numreads"] / (
        df_amr["Reference.sequence.length"] * df_amr["total_reads"]
    )
    df_amr["rpkm2"] = 1e9 * df_amr["numreads"] / (
        df_amr["Alignment.length"] * df_amr["total_reads"]
    )
    
    # FPKM using mapped fraction as per your R script
    df_amr["fpkm"] = 1e9 * df_amr["numreads"] / (
        df_amr["Reference.sequence.length"]
        * (df_amr["total_reads"] * (df_amr["overall_aln_percent"] / 100.0))
    )
    
    # FPKM_prok only if SingleM SPF is being used
    if spf_files_pe:
        # extract prokaryotic base nums combine across pe and se singlem outputs
        prok_bases_pe = pd.to_numeric(df_amr["bacterial_archaeal_bases_pe"], errors="coerce").fillna(0)
        prok_bases_s1 = pd.to_numeric(df_amr["bacterial_archaeal_bases_s1"], errors="coerce").fillna(0)
        prok_bases_s2 = pd.to_numeric(df_amr["bacterial_archaeal_bases_s2"], errors="coerce").fillna(0)
        prok_bases_tot = prok_bases_pe + prok_bases_s1 + prok_bases_s2
        
        total_bases_pe = pd.to_numeric(df_amr["metagenome_size_pe"], errors="coerce").fillna(0)
        total_bases_s1 = pd.to_numeric(df_amr["metagenome_size_s2"], errors="coerce").fillna(0)
        total_bases_s2 = pd.to_numeric(df_amr["metagenome_size_s2"], errors="coerce").fillna(0)
        total_bases_tot = total_bases_pe + total_bases_s1 + total_bases_s2
        
        df_amr["Reference.sequence.length"] = pd.to_numeric(df_amr["Reference.sequence.length"], errors="coerce")
        
        # obtain SPF as a proportion (read_fraction), not percentage 
        prok_frac = prok_bases_tot / total_bases_tot
        eff_lib = df_amr["total_reads"] * prok_frac
        df_amr["fpkm_prok"] = (1e9 * df_amr["numreads"]) / (
            df_amr["Reference.sequence.length"] * eff_lib
        )
        df_amr.loc[(eff_lib <= 0) | eff_lib.isna(), "fpkm_prok"] = pd.NA
    else:
        df_amr["fpkm_prok"] = pd.NA
    
    # Write output
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_amr.to_csv(out_csv, index=False)
    eprint(f"Wrote: {out_csv}")
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
