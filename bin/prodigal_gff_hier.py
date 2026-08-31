#!/usr/bin/env python3
"""Re-emit a Prodigal CDS-only GFF3 as full gene->mRNA->CDS blocks.

Prodigal's native `-f gff` output contains CDS features only (no gene/mRNA
rows), but EVM assembles consensus models from gene/mRNA blocks: a CDS-only
pass-through source is structurally inert no matter its weight. On the
microsporidia OC4 validation genome the CDS-only Prodigal track was silently
ignored and the consensus degraded to a pure GeneMark echo (gene Sn 0.584);
with this transform feeding funannotate predict's --other_gff, gene Sn/Sp
recovered to 0.839/0.817 vs 0.849/0.791 for Prodigal alone.

Each CDS becomes a gene/mRNA/CDS block separated by a blank line (this is
what funannotate-runEVM.py's lib.readBlocks / predict.py's blank-insertion
normalization expect). No ##gff-version header is written: a header line
joins the first block and crashes EVM's gene_blocks_to_interlap, which
expects every block to be immediately preceded by a blank line.

Usage:
    prodigal_gff_hier.py PRODIGAL.gff3 -o OUT.gff3
"""
import argparse
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gff3", help="Prodigal -f gff output (CDS-only)")
    parser.add_argument("-o", "--output", required=True, help="Output GFF3 path")
    args = parser.parse_args()

    n = 0
    with open(args.gff3) as fh, open(args.output, "w") as o:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[2] != "CDS":
                continue
            seqid, _, _, start, end, score, strand, phase, _ = cols
            n += 1
            gid, mid = f"prodigal_g{n}", f"prodigal_m{n}"
            o.write(f"{seqid}\tprodigal\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gid}\n")
            o.write(
                f"{seqid}\tprodigal\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={mid};Parent={gid}\n"
            )
            o.write(
                f"{seqid}\tprodigal\tCDS\t{start}\t{end}\t{score}\t{strand}\t{phase}\tID={mid}.cds;Parent={mid}\n"
            )
            o.write("\n")

    if n == 0:
        print(f"ERROR: no CDS features found in {args.gff3}", file=sys.stderr)
        return 1
    print(f"wrote {n} gene/mRNA/CDS blocks to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
