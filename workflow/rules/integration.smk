configfile: "config/integration.yaml"

rule integration_atac_benchmark_features:
    input:
        benchmark_tsv=config["integration"]["benchmark"]["gold_tsv"],
        benchmark_bed=config["integration"]["benchmark"]["gold_bed"],
        idr10=config["integration"]["atac"]["idr10"],
        idr05=config["integration"]["atac"]["idr05"]

    output:
        features="results/integration/atac/benchmark_atac_features.tsv",
        anchors="results/integration/atac/benchmark_gold_anchors.bed",
        idr10_unique="results/integration/atac/atac_idr10_unique.bed",
        idr05_unique="results/integration/atac/atac_idr05_unique.bed",
        summary="results/integration/atac/integration_summary.tsv",
        provenance="results/integration/atac/benchmark_atac_features.provenance.json"

    params:
        windows=lambda wildcards: ",".join(
            str(x)
            for x in config["integration"]["atac"]["anchor_half_windows_bp"]
        )

    threads: 1

    shell:
        r"""
        python workflow/scripts/build_atac_benchmark_features.py \
          --benchmark-tsv {input.benchmark_tsv:q} \
          --benchmark-bed {input.benchmark_bed:q} \
          --idr10 {input.idr10:q} \
          --idr05 {input.idr05:q} \
          --windows {params.windows:q} \
          --features {output.features:q} \
          --anchors-bed {output.anchors:q} \
          --idr10-unique-bed {output.idr10_unique:q} \
          --idr05-unique-bed {output.idr05_unique:q} \
          --summary {output.summary:q} \
          --provenance {output.provenance:q}
        """
