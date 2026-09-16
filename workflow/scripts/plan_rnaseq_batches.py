#!/usr/bin/env python3
import argparse, csv, statistics, json
from selection_gate import ROOT, decision, catalogue, fingerprint
from collections import defaultdict
from pathlib import Path

def read(path):
    with path.open(encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))

def write(path, rows, columns=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]) if rows else columns, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

def fit(xs, ys):
    mx, my = statistics.mean(xs), statistics.mean(ys)
    b = sum((x-mx)*(y-my) for x,y in zip(xs,ys))/sum((x-mx)**2 for x in xs)
    return my-b*mx, b

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--benchmark-metrics", type=Path, required=True)
    p.add_argument("--run-sizes", type=Path, default=Path("results/rnaseq/benchmark/rnaseq_run_sizes.tsv"))
    p.add_argument("--samples", type=Path, default=Path("config/samples.tsv"))
    p.add_argument("--classes", type=Path, default=Path("config/rnaseq_resource_classes.tsv"))
    p.add_argument("--outdir", type=Path, default=Path("results/rnaseq/benchmark"))
    a = p.parse_args()
    metrics, runs, class_rows = read(a.benchmark_metrics), read(a.run_sizes), read(a.classes)
    if a.samples.resolve() != (ROOT / "config/samples.tsv").resolve():
        raise SystemExit("--samples must be this checkout's config/samples.tsv for shared selection provenance")
    platform_by_run = {acc: row.get("instrument_platform", "").strip().upper() for acc,row in catalogue().items()}
    xs = [int(r["fastq_bytes"])/1e9 for r in metrics]
    ti, ts = fit(xs, [int(r["elapsed_seconds"])/60 for r in metrics])
    si, ss = fit(xs, [int(r["incremental_peak_tmpdir_bytes"])/1e9 for r in metrics])
    baseline = statistics.median(int(r["baseline_tmpdir_bytes"])/1e9 for r in metrics)
    classes = []
    for i, r in enumerate(class_rows):
        classes.append({**r, "min": int(r["min_fastq_bytes_exclusive"]),
                        "max": int(r["max_fastq_bytes_inclusive"]) if r["max_fastq_bytes_inclusive"] else None,
                        "capacity": float(r["batch_target_predicted_minutes"])})
    def classify(size):
        for i,c in enumerate(classes):
            if (size >= c["min"] if i == 0 else size > c["min"]) and (c["max"] is None or size <= c["max"]):
                return c
        raise ValueError(f"No class for {size}")
    plan, grouped = [], defaultdict(list)
    for r in runs:
        size = int(r["fastq_bytes"]); gb = size/1e9; c = classify(size)
        structure = (r.get("fastq_structure") or r.get("processing_layout") or r.get("library_layout") or "UNKNOWN").strip()
        platform = platform_by_run.get(r["run_accession"], "")
        eligible, selection_reason = decision(r["run_accession"], omics="RNA-seq")
        scheduling_status = (
            "SCHEDULED"
            if eligible
            else "BLOCKED_SELECTION"
        )
        item = {"run_accession": r["run_accession"], "study_accession": r.get("study_accession", ""),
                "fastq_structure": structure, "instrument_platform": platform, "fastq_bytes": size, "fastq_GB": f"{gb:.6f}",
                "resource_class": c["resource_class"], "predicted_minutes": f"{max(1,ti+ts*gb):.3f}",
                "predicted_scratch_apparent_GB": f"{baseline+max(0,si+ss*gb):.3f}",
                "cpus": c["cpus"], "memory_mb": c["memory_mb"], "scratch_target_gb": c["scratch_target_gb"],
                "estimated_allocated_tmpdir_gb": f"{845*int(c['cpus'])/64:.3f}",
                "slurm_constraint": c["slurm_constraint"], "walltime": c["walltime"],
                "evidence_level": c["evidence_level"],
                "scheduling_status": scheduling_status, "selection_reason": selection_reason,
                "batch_id": "", "batch_order": ""}
        plan.append(item)
        if item["scheduling_status"] == "SCHEDULED": grouped[c["resource_class"]].append(item)
    batch_rows, summaries = [], []
    for c in classes:
        name, bins = c["resource_class"], []
        for item in sorted(grouped[name], key=lambda x: (-float(x["predicted_minutes"]), x["run_accession"])):
            candidates = [(c["capacity"]-(b["minutes"]+float(item["predicted_minutes"])), i)
                          for i,b in enumerate(bins) if b["minutes"]+float(item["predicted_minutes"]) <= c["capacity"]]
            if candidates:
                _, i = min(candidates); bins[i]["runs"].append(item); bins[i]["minutes"] += float(item["predicted_minutes"])
            else:
                bins.append({"minutes": float(item["predicted_minutes"]), "runs": [item]})
        for n,b in enumerate(bins,1):
            bid = f"{name}_{n:03d}"
            ordered = sorted(b["runs"], key=lambda x: (-float(x["predicted_minutes"]), x["run_accession"]))
            for order,item in enumerate(ordered,1):
                item["batch_id"], item["batch_order"] = bid, order
                batch_rows.append({"batch_id": bid, "batch_order": order, "resource_class": name,
                                   "run_accession": item["run_accession"], "study_accession": item["study_accession"],
                                   "fastq_structure": item["fastq_structure"], "instrument_platform": item["instrument_platform"], "fastq_bytes": item["fastq_bytes"],
                                   "predicted_minutes": item["predicted_minutes"]})
            summaries.append({"batch_id": bid, "resource_class": name, "n_runs": len(ordered),
                              "total_fastq_GB": f"{sum(int(x['fastq_bytes']) for x in ordered)/1e9:.3f}",
                              "predicted_processing_minutes": f"{b['minutes']:.3f}",
                              "maximum_predicted_scratch_apparent_GB": f"{max(float(x['predicted_scratch_apparent_GB']) for x in ordered):.3f}",
                              "cpus": c["cpus"], "memory_mb": c["memory_mb"], "scratch_target_gb": c["scratch_target_gb"],
                              "estimated_allocated_tmpdir_gb": f"{845*int(c['cpus'])/64:.3f}",
                              "slurm_constraint": c["slurm_constraint"], "walltime": c["walltime"],
                              "evidence_level": c["evidence_level"]})
    write(a.outdir/"rnaseq_execution_plan.tsv", sorted(plan, key=lambda x:(x["scheduling_status"],x["resource_class"],x["batch_id"],int(x["batch_order"] or 0))))
    write(a.outdir/"rnaseq_batches.tsv", batch_rows, ["batch_id","batch_order","resource_class","run_accession","study_accession","fastq_structure","instrument_platform","fastq_bytes","predicted_minutes"])
    write(a.outdir/"rnaseq_batch_summary.tsv", summaries, ["batch_id","resource_class","n_runs"])
    (a.outdir/"selection_provenance.json").write_text(json.dumps(fingerprint(), indent=2)+"\n")
    scheduled = [x for x in plan if x["scheduling_status"] == "SCHEDULED"]
    if len(plan) != len({x["run_accession"] for x in plan}) or len(batch_rows) != len(scheduled):
        raise SystemExit("ERROR: duplicated or unassigned runs")
    print(f"Total runs represented: {len(plan)}\nScheduled supported short-read runs: {len(scheduled)}\nBlocked selection/platform runs: {len(plan)-len(scheduled)}\nTotal batches: {len(summaries)}")
    print("\nclass\tscheduled\texcluded\tbatches")
    for c in classes:
        n=c["resource_class"]
        print(n, sum(x["resource_class"]==n for x in scheduled), sum(x["resource_class"]==n and x["scheduling_status"]!="SCHEDULED" for x in plan), sum(x["resource_class"]==n for x in summaries), sep="\t")
    print(f"\nTime model: minutes = {ti:.6f} + {ts:.6f} x FASTQ_GB")
    print(f"Scratch model: apparent_GB = {baseline+si:.6f} + {ss:.6f} x FASTQ_GB")

if __name__ == "__main__": main()
