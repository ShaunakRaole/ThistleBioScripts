#!/usr/bin/env python3
import sys
import argparse
import pandas as pd
import re

def parse_attribute(attr_str):
    """
    Parse a GTF attribute string into a dict,
    matching both gene_id "XYZ";  and  gene_id XYZ;
    """
    d = {}
    # This will match either "XYZ" or XYZ (up to the semicolon)
    for m in re.finditer(r'(\S+)\s+(?:"([^"]+)"|([^";\s]+));', attr_str):
        key = m.group(1)
        val = m.group(2) if m.group(2) is not None else m.group(3)
        d[key] = val
    return d

def main(in_gtf, out_gtf):
    # --- 1) Read GTF and subset to exons ---
    cols = ["seqname","source","feature","start","end","score","strand","frame","attribute"]
    df = pd.read_csv(
        in_gtf, sep="\t", comment="#", header=None, names=cols,
        dtype={"start":int,"end":int}, keep_default_na=False
    )
    ex = df[df.feature == "exon"].copy()
    
    # Parse out gene_id & transcript_id from column 9
    attrs = ex.attribute.apply(parse_attribute)
    attrs = ex.attribute.apply(parse_attribute)
    ex["gene_id"] = attrs.apply(lambda d: d.get("gene_id",""))
    ex["transcript_id"] = attrs.apply(lambda d: d.get("transcript_id",""))
    # fill any empty transcript_id with the gene_id
    ex["transcript_id"] = ex.apply(
        lambda r: r.gene_id if not r.transcript_id else r.transcript_id,
        axis=1
    )
    
    # --- 2) Sort so that we can sweep merge/nest-drop by chr+strand, start asc ---
    ex = ex.sort_values(
        ["seqname","strand","start","end"],
        ascending=[True,True,True,False]
    ).reset_index(drop=True)
    
    # --- 3) Sweep through and drop fully nested exons, merge partial overlaps ---
    merged_blocks = []
    current = None

    for _, row in ex.iterrows():
        s, e = row.start, row.end
        gid, tid = row.gene_id, row.transcript_id
        src = row.source
        strand = row.strand
        seq = row.seqname
        
        if current is None:
            # start first block
            current = {
                "seqname": seq,
                "source":  src,
                "strand":  strand,
                "start":   s,
                "end":     e,
                "max_len": e - s + 1,
                "gene_id": gid,
                "transcript_id": tid
            }
        else:
            # same chr+strand?
            if seq == current["seqname"] and strand == current["strand"] and s <= current["end"]:
                # overlap or nested → extend end if needed
                if e > current["end"]:
                    current["end"] = e
                # pick the exon with the longer original span
                length = e - s + 1
                if length > current["max_len"]:
                    current["max_len"]       = length
                    current["gene_id"]       = gid
                    current["transcript_id"] = tid
            else:
                # disjoint → flush current block
                merged_blocks.append(current)
                # start new
                current = {
                    "seqname": seq,
                    "source":  src,
                    "strand":  strand,
                    "start":   s,
                    "end":     e,
                    "max_len": e - s + 1,
                    "gene_id": gid,
                    "transcript_id": tid
                }
    # flush last
    if current is not None:
        merged_blocks.append(current)
    
    # --- 4) Write out merged exons as a GTF ---
    with open(out_gtf, "w") as out:
        for b in merged_blocks:
            # build attribute string
            attr = f'gene_id "{b["gene_id"]}"; transcript_id "{b["transcript_id"]}";'
            cols9 = [b["seqname"], b["source"], "exon",
                     str(b["start"]), str(b["end"]),
                     ".", b["strand"], ".", attr]
            out.write("\t".join(cols9) + "\n")

if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Drop nested exons and merge partial overlaps (keep longest gene_id)."
    )
    p.add_argument("input_gtf",  help="Sorted GTF with individual exon entries")
    p.add_argument("output_gtf", help="Output GTF of pruned+merged exons")
    args = p.parse_args()
    main(args.input_gtf, args.output_gtf)