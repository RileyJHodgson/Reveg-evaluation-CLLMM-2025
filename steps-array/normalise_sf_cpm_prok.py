#!/usr/bin/env python3

import json
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(
        description="CPM normalisation of functional table using prokaryotic fraction"
    )
    parser.add_argument("--qc-json", required=True, help="fastp JSON file")
    parser.add_argument("--spf-file", required=True, help="singlem SPF TSV file")
    parser.add_argument("--sf-file", required=True, help="input functional table")
    parser.add_argument("--out-file", required=True, help="output CPM-normalised table")

    args = parser.parse_args()

    # --- extract total reads ---
    with open(args.qc_json) as f:
        qc = json.load(f)

    total_reads = qc["summary"]["after_filtering"]["total_reads"]

    # --- extract prokaryotic fraction ---
    with open(args.spf_file) as f:
        lines = f.readlines()

    if len(lines) < 2:
        raise RuntimeError("SPF file does not contain a second line")

    prok_frac = float(lines[1].rstrip().split("\t")[3])

    print(total_reads)
    print(prok_frac)

    # --- denominator ---
    denom = total_reads * (prok_frac / 100.0)

    if denom <= 0:
        raise RuntimeError("Computed denominator is zero or negative")

    # --- apply CPM normalisation ---
    with open(args.sf_file) as fin, open(args.out_file, "w") as fout:
        for line_num, line in enumerate(fin, start=1):
            line = line.rstrip("\n")

            # Metadata lines (1–4)
            if line_num <= 4:
                fout.write(line + "\n")
                continue

            # Header line (5)
            if line_num == 5:
                fout.write(line + "\n")
                continue

            # Data lines (6+)
            fields = line.split("\t")

            for i in range(4, len(fields)):  # columns 5+ (0-based index)
                if fields[i] not in ("", "NA"):
                    try:
                        value = float(fields[i])
                        fields[i] = f"{(value / denom) * 1e6:.6f}"
                    except ValueError:
                        # Leave non-numeric fields unchanged
                        pass

            fout.write("\t".join(fields) + "\n")

    print("CPM_prok normalisation complete")
    print(f"Output: {args.out_file}")


if __name__ == "__main__":
    main()
