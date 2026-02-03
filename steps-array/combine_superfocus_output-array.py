#!/usr/bin/env python3

import argparse
import glob
import os
import re
import sys
from functools import reduce

import pandas as pd

KEYS = ["subsystem_1", "subsystem_2", "subsystem_3", "functional_process"]

def read_superfocus_table(path: str) -> pd.DataFrame:
    """
    Read a Superfocus *_all_levels_and_function.xls file that is actually a
    tab-delimited text table with 4 metadata lines before the header.
    Mirrors: read_tsv(file, skip = 4).
    """
    # skip first 4 lines; line 5 is header
    df = pd.read_csv(path, sep="\t", skiprows=4, dtype=str)
    
    # After read, last two columns should be counts and perc (but their names can vary)
    if df.shape[1] < 6:
        raise ValueError(f"Expected >=6 columns after parsing {path}, got {df.shape[1]}")
    
    return df


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Combine and reformat Superfocus *_all_levels_and_function.xls outputs."
    )
    parser.add_argument(
        "--input-dir",
        default=".",
        help="Directory to search (recursively) for *_all_levels_and_function.xls files (default: .)",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Where to write outputs. Default: <input-dir>/OTU_style_output and <input-dir>/Long_format_output",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search recursively under input-dir (recommended for pipeline outputs).",
    )
    args = parser.parse_args()
    
    input_dir = os.path.abspath(args.input_dir)
    
    if args.output_dir is None:
        main_dir = input_dir
    else:
        main_dir = os.path.abspath(args.output_dir)
    
    out_dir = main_dir
    
    pattern = "**/func/*_all_levels_and_function.xls" if args.recursive else "*_all_levels_and_function.xls"
    files = sorted(glob.glob(os.path.join(input_dir, pattern), recursive=args.recursive))
    
    if not files:
        print(f"ERROR: No files matched pattern {pattern} under {input_dir}", file=sys.stderr)
        return 1
    
    sample_dfs = []
    suffix_re = re.compile(r"_all_levels_and_function\.xls$")
    
    for f in files:
        base = os.path.basename(f)
        if not suffix_re.search(base):
            continue
        
        sample = suffix_re.sub("", base)
        
        df = read_superfocus_table(f)
        
        # Rename first 4 columns to standard keys. Rename last 2 columns to sample and sample_PERC.
        # If your file has extra columns beyond 6, we keep them but you may want to drop them explicitly.
        cols = list(df.columns)
        
        # Standardise key columns
        cols[0:4] = KEYS
        
        # Standardise the final two columns as sample count and sample perc
        cols[-2] = sample
        cols[-1] = f"{sample}_PERC"
        
        df.columns = cols
        
        # Keep just keys + last two columns (mirrors typical Superfocus table structure)
        df = df[KEYS + [sample, f"{sample}_PERC"]]
        
        # Convert numeric columns; non-numeric -> NaN -> will become 0 later
        df[sample] = pd.to_numeric(df[sample], errors="coerce")
        df[f"{sample}_PERC"] = pd.to_numeric(df[f"{sample}_PERC"], errors="coerce")
        
        sample_dfs.append(df)
        print(f"Loaded: {f}")
        
    if not sample_dfs:
        print("ERROR: No valid Superfocus files were read.", file=sys.stderr)
        return 1
        
    # Outer merge across all samples on KEYS (equivalent to Reduce(merge(..., all=TRUE)) )
    merged = reduce(lambda left, right: pd.merge(left, right, on=KEYS, how="outer"), sample_dfs)
    
    # Replace NA with 0 for numeric columns; keys remain strings
    for c in merged.columns:
        if c not in KEYS:
            merged[c] = merged[c].fillna(0)

    # Count-only (drop *_PERC)
    count_cols = [c for c in merged.columns if not c.endswith("_PERC")]
    merged_count = merged[count_cols]
    
    # Perc-only (keys + *_PERC)
    perc_cols = KEYS + [c for c in merged.columns if c.endswith("_PERC")]
    merged_perc = merged[perc_cols]
    
    merged_count.to_csv(os.path.join(out_dir, "functions_counts.csv"), index=False)
    merged_perc.to_csv(os.path.join(out_dir, "functions_percentages.csv"), index=False)

    print(f"Wrote:\n  {os.path.join(out_dir, 'functions_counts.csv')}\n"
          f"  {os.path.join(out_dir, 'functions_percentages.csv')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

