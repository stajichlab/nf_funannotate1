#!/usr/bin/env python3
"""Merge per-task *.versions.yml files (from record_versions.sh) into TSVs.

Usage: merge_versions.py ROOT [--all OUT.tsv]

ROOT holds <sample>/versions/*.versions.yml. Writes
ROOT/<sample>/versions/<sample>.tool_versions.tsv per sample, and optionally
one combined TSV. Columns: sample, step, tool, version, container, pipeline,
host, recorded, retroactive, note.
"""
import argparse
import glob
import os
import re

COLS = ["sample", "step", "tool", "version", "container", "pipeline", "host", "recorded", "retroactive", "note"]
KV = re.compile(r'^(\s*)([^:\s][^:]*):\s*"(.*)"\s*$')


def parse(path):
    meta, tools, proc, in_tools = {}, [], "", False
    with open(path) as fh:
        for line in fh:
            if line.startswith('"'):
                proc = line.strip().rstrip(":").strip('"')
                continue
            m = KV.match(line)
            if not m:
                in_tools = in_tools or line.strip() == "tools:"
                continue
            val = m.group(3).replace('\\"', '"').replace("\\\\", "\\")
            if len(m.group(1)) >= 4:
                tools.append((m.group(2), val))
            else:
                meta[m.group(2)] = val
    return proc, meta, tools


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--all", dest="all_out")
    a = ap.parse_args()
    rows_all = []
    for d in sorted(glob.glob(os.path.join(a.root, "*", "versions"))):
        sample = os.path.basename(os.path.dirname(d))
        rows = []
        for f in sorted(glob.glob(os.path.join(d, "*.versions.yml"))):
            proc, meta, tools = parse(f)
            step = proc.split(":")[-1]
            for tool, ver in tools:
                rows.append([sample, step, tool, ver] + [meta.get(k, "") for k in COLS[4:]])
        with open(os.path.join(d, f"{sample}.tool_versions.tsv"), "w") as out:
            out.write("\t".join(COLS) + "\n")
            out.writelines("\t".join(r) + "\n" for r in rows)
        rows_all += rows
    if a.all_out:
        with open(a.all_out, "w") as out:
            out.write("\t".join(COLS) + "\n")
            out.writelines("\t".join(r) + "\n" for r in rows_all)
    print(f"{len(rows_all)} records")


if __name__ == "__main__":
    main()
