#!/usr/bin/env python3
"""Technical validation of one eligible SINGLE ATAC production result."""
import argparse,csv,json,math
from pathlib import Path
from atac_manifest import eligible_run
from file_hash import sha256_file


def table(path):
    with open(path,newline="",encoding="utf-8") as handle:return list(csv.DictReader(handle,delimiter="\t"))


def validate_tss(profile,summary,window,bin_size,background):
    rows=table(summary)
    if not rows or any(set(r)!={"metric","value"} for r in rows):raise ValueError("Malformed TSS summary")
    if len({r["metric"] for r in rows})!=len(rows):raise ValueError("Duplicate TSS metric")
    try:values={r["metric"]:float(r["value"]) for r in rows}
    except ValueError as exc:raise ValueError("Nonnumeric TSS summary") from exc
    required={"usable_tss","window_bp","bin_size_bp","background_bp_per_side","background_mean","tss_enrichment_score"}
    if required-set(values) or any(not math.isfinite(v) or v<0 for v in values.values()):raise ValueError("Invalid TSS summary values")
    if (values["window_bp"],values["bin_size_bp"],values["background_bp_per_side"])!=(window,bin_size,background):raise ValueError("TSS parameters disagree")
    bins=table(profile)
    if len(bins)!=2*window//bin_size:raise ValueError("Unexpected TSS profile length")
    for i,row in enumerate(bins):
        try:start=int(row["relative_start"]);end=int(row["relative_end"]);count=int(row["insertion_count"]);signal=float(row["normalized_signal"])
        except (KeyError,ValueError) as exc:raise ValueError("Malformed TSS profile") from exc
        if start!=-window+i*bin_size or end!=start+bin_size or count<0 or not math.isfinite(signal) or signal<0:raise ValueError("Invalid TSS profile bin")
    return values


def main():
    import pysam
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("run","manifest","bam","bai","profile","summary","peaks","frip","bigwig","output"):p.add_argument("--"+name,required=True)
    p.add_argument("--min-mapq",type=int,default=30);p.add_argument("--mitochondrial",default="NC_007936.1")
    p.add_argument("--window",type=int,default=2000);p.add_argument("--bin-size",type=int,default=10);p.add_argument("--background",type=int,default=100)
    a=p.parse_args();eligible_run(a.manifest,a.run,require_single=True)
    count=0;pysam.quickcheck(a.bam)
    with pysam.AlignmentFile(a.bam,"rb",index_filename=a.bai,require_index=True) as bam:
        if bam.header.to_dict().get("HD",{}).get("SO")!="coordinate":raise ValueError("Filtered BAM is not coordinate sorted")
        for read in bam.fetch(until_eof=True):
            count+=1
            if read.flag&3844 or read.mapping_quality<a.min_mapq or read.reference_name==a.mitochondrial or read.is_paired:
                raise ValueError("Filtered BAM violates SINGLE filtering contract")
    if count==0:raise ValueError("Empty filtered BAM")
    tss=validate_tss(a.profile,a.summary,a.window,a.bin_size,a.background)
    paths={k:getattr(a,k) for k in ("bam","bai","profile","summary","peaks","frip","bigwig")}
    value={"run_accession":a.run,"technical_validation":"passed","biological_acceptance":"not_assessed",
      "filtered_records":count,"tss_metrics":tss,"sha256":{k:sha256_file(v) for k,v in paths.items()}}
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(value,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()
