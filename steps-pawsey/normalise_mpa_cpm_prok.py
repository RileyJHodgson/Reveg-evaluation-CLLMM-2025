#!/usr/bin/env python3

import json
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(
        description="Normalize mpa breport to CPM using prokaryotic fraction"
    )
    parser.add_argument("--qc-json", required=True, help="fastp JSON file")
    parser.add_argument("--spf-file", required=True, help="singlem spf TSV")
    parser.add_argument("--mpa-file", required=True, help="mpa breport input")
    parser.add_argument("--out-file", required=True, help="output TSV")

    args = parser.parse_args()

    # --- extract total reads from fastp JSON ---
    with open(args.qc_json) as f:
        qc = json.load(f)

    total_reads = qc["summary"]["after_filtering"]["total_reads"]
    print(total_reads)
    # --- extract prokaryotic fraction ---
    with open(args.spf_file) as f:
        lines = f.readlines()

    if len(lines) < 2:
        raise RuntimeError("SPF file does not contain expected data")

    prok_frac = float(lines[1].rstrip().split("\t")[3])
    print(prok_frac)

    # --- denominator ---
    denom = total_reads * (prok_frac / 100.0)
    print(denom)

    if denom <= 0:
        raise RuntimeError("Denominator is zero or negative")

    # --- normalize mpa file ---
    with open(args.mpa_file) as fin, open(args.out_file, "w") as fout:
        for i, line in enumerate(fin, start=1):
            line = line.rstrip("\n")

            # header
            if i == 1:
                fout.write(line + "\tCPM_prok\n")
                continue

            fields = line.split("\t")

            try:
                value = float(fields[1])
                cpm = (value / denom) * 1e6
                fields[1] = f"{cpm:.6f}"
            except (IndexError, ValueError):
                # leave line unchanged if malformed
                pass

            fout.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    main()
