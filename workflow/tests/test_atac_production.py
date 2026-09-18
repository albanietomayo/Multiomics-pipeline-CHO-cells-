"""Small deterministic ATAC production contracts; no biological validation."""
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"workflow/scripts"))
from atac_manifest import eligible_run,eligible_runs,load_groups
from atac_reference_span import calculate as calculate_reference_span
from build_atac_production_manifest import build_rows,write_manifest
from calculate_atac_frip import count_frip
from calculate_atac_tss_enrichment import load_tss_bed,usable_read
from clip_atac_bedgraph import clip
from filter_atac_bam import filter_bam,rejected,validate_pair,nuclear_references
from filter_idr import filter_idr
from file_hash import sha256_file
from tn5_transform import insertion_position,transform
from validate_atac_bigwig import signal_chromosomes,validate_bigwig,validate_header
from validate_atac_run import validate_tss


class Read:
    def __init__(self,name="x",flag=0,mapq=30,chrom="nuclear",rid=0,start=10,mrid=0,mstart=20,read1=False,read2=False):
        self.query_name=name;self.flag=flag;self.mapping_quality=mapq;self.reference_name=chrom
        self.reference_id=rid;self.reference_start=start;self.next_reference_id=mrid;self.next_reference_start=mstart
        self.is_read1=read1;self.is_read2=read2


class ATACProductionContracts(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.d=Path(self.temp.name)
        self.manifest=self.d/"atac.tsv"
        self.source_manifest=ROOT/"workflow/tests/fixtures/atac_fastq_manifest.tsv"
        rows=build_rows(self.source_manifest,ROOT)
        write_manifest(rows,self.manifest)

    def tearDown(self):self.temp.cleanup()

    def test_production_eligibility_exact_regression(self):
        self.assertEqual(eligible_runs(self.manifest),["SRR12774931","SRR12774932"])
        for run in ("SRR12774931","SRR12774932"):eligible_run(self.manifest,run)
        for run in ("SRR12774934","SRR29929613","SRR29929615","SRR29929617"):
            with self.subTest(run=run),self.assertRaises(ValueError):eligible_run(self.manifest,run)

    def test_paired_downstream_fails_closed_even_if_selection_eligible(self):
        rows=[]
        with self.manifest.open() as handle:
            for row in csv.DictReader(handle,delimiter="\t"):
                if row["run_accession"]=="SRR29929613":
                    row["selection_eligible"]="true";row["selection_decision"]="eligible";row["selection_reason"]="eligible"
                    rows.append(row)
        paired=self.d/"paired_eligible.tsv"
        with paired.open("w",newline="") as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n")
            writer.writeheader();writer.writerows(rows)
        eligible_run(paired,"SRR29929613")
        with self.assertRaisesRegex(ValueError,"validated only for SINGLE"):
            eligible_run(paired,"SRR29929613",require_single=True)

    def test_manifest_preserves_file_level_fields_and_fails_closed(self):
        with self.manifest.open() as handle:
            rows=list(csv.DictReader(handle,delimiter="\t"))
        self.assertEqual(len(rows),9)
        required={"study_accession","sample_accession","library_layout","processing_layout","fastq_structure","fastq_role","selection_decision","selection_reason"}
        self.assertFalse(required-set(rows[0]))
        source=self.source_manifest
        bad=self.d/"bad.tsv"
        with source.open() as handle:
            data=list(csv.DictReader(handle,delimiter="\t"))
        fields=list(data[0])
        for row in data:
            if row["run_accession"]=="SRR12774931":row["fastq_role"]="UNRESOLVED"
        with bad.open("w",newline="") as h:w=csv.DictWriter(h,fields,delimiter="\t");w.writeheader();w.writerows(data)
        built=build_rows(bad,ROOT);q=[r for r in built if r["run_accession"]=="SRR12774931"]
        self.assertEqual(q[0]["selection_eligible"],"false");self.assertIn("inconsistent_fastq_roles",q[0]["selection_reason"])

    def test_filter_flags_mapq_mito_and_pair_membership(self):
        self.assertFalse(rejected(Read(),30,3844,"NC_007936.1"))
        for read in (Read(flag=1024),Read(flag=256),Read(flag=2048),Read(flag=512),Read(mapq=29),Read(chrom="NC_007936.1")):
            with self.subTest(flag=read.flag,mapq=read.mapping_quality,chrom=read.reference_name):
                self.assertTrue(rejected(read,30,3844,"NC_007936.1"))
        a=Read("p",flag=99,start=10,mstart=20,read1=True);b=Read("p",flag=147,start=20,mstart=10,read2=True)
        self.assertTrue(validate_pair([a,b]))
        self.assertFalse(validate_pair([a]));self.assertFalse(validate_pair([a,Read("q",flag=147,start=20,mstart=10,read2=True)]))
        b.next_reference_start=999;self.assertFalse(validate_pair([a,b]))

    def make_bam(self,records,references=(("nuclear",1000),("NC_007936.1",500))):
        import pysam
        path=self.d/("input_"+str(len(list(self.d.glob("input_*.bam"))))+".bam")
        header={"HD":{"VN":"1.6","SO":"coordinate"},"SQ":[{"SN":n,"LN":s} for n,s in references]}
        with pysam.AlignmentFile(path,"wb",header=header) as bam:
            for name,flag,rid,start,mapq,mrid,mstart in records:
                read=pysam.AlignedSegment();read.query_name=name;read.query_sequence="A"*30
                read.flag=flag;read.reference_id=rid;read.reference_start=start;read.mapping_quality=mapq;read.cigarstring="30M"
                if flag&1:read.next_reference_id=mrid;read.next_reference_start=mstart
                bam.write(read)
        return path

    def test_real_synthetic_single_filter(self):
        import pysam
        source=self.make_bam([("good",0,0,10,30,-1,-1),("low",0,0,50,29,-1,-1),
          ("secondary",256,0,100,40,-1,-1),("supplementary",2048,0,150,40,-1,-1),
          ("qcfail",512,0,200,40,-1,-1),("duplicate",1024,0,250,40,-1,-1),("mito",0,1,10,40,-1,-1)])
        output=self.d/"single.bam";metrics=self.d/"single.json"
        self.assertEqual(filter_bam(source,output,"SINGLE",30,3844,3,3852,"NC_007936.1",self.d,metrics),1)
        with pysam.AlignmentFile(output,"rb") as bam:self.assertEqual([r.query_name for r in bam.fetch(until_eof=True)],["good"])
        self.assertEqual(json.loads(metrics.read_text())["output_records"],1)

    def test_real_synthetic_paired_filter_and_orphans(self):
        import pysam
        records=[("good",99,0,10,40,0,100),("orphan",99,0,50,40,0,150),
          ("good",147,0,100,40,0,10),("orphan",147,0,150,10,0,50)]
        source=self.make_bam(records);output=self.d/"paired.bam"
        self.assertEqual(filter_bam(source,output,"PAIRED",30,3844,3,3852,"NC_007936.1",self.d,self.d/"paired.json"),2)
        with pysam.AlignmentFile(output,"rb") as bam:self.assertEqual([r.query_name for r in bam.fetch(until_eof=True)],["good","good"])
        bad=self.make_bam([("a",99,0,10,40,0,100),("b",147,0,100,40,0,10)])
        with self.assertRaises(ValueError):filter_bam(bad,self.d/"bad.bam","PAIRED",30,3844,3,3852,"NC_007936.1",self.d,self.d/"bad.json")

    def test_real_tn5_representations(self):
        source=self.make_bam([("plus",0,0,10,40,-1,-1),("minus",16,0,100,40,-1,-1)])
        insertions=self.d/"insertions.bed";tags=self.d/"tags.bed";metrics=self.d/"tn5.json"
        value=transform(source,insertions,tags,metrics)
        self.assertEqual(insertions.read_text().splitlines(),[
          "nuclear\t14\t15\tplus\t40\t+","nuclear\t125\t126\tminus\t40\t-"])
        self.assertEqual(tags.read_text().splitlines(),[
          "nuclear\t14\t40\tplus\t40\t+","nuclear\t100\t125\tminus\t40\t-"])
        self.assertEqual((value["input_records"],value["shifted_intervals"]),(2,2))

    def test_tn5_transform_externally_sorts_shifted_intervals(self):
        # BAM order is plus@9 then reverse@10, but +4 makes their tagAlign
        # starts 13 then 10. The transform must externally re-sort that output.
        source=self.make_bam([("plus",0,0,9,40,-1,-1),("minus",16,0,10,40,-1,-1)])
        insertions=self.d/"ordered.insertions.bed";tags=self.d/"ordered.tags.bed"
        transform(source,insertions,tags,self.d/"ordered.metrics.json",tmpdir=self.d)
        self.assertEqual(tags.read_text().splitlines(),[
          "nuclear\t10\t35\tminus\t40\t-","nuclear\t13\t39\tplus\t40\t+"])
        self.assertEqual(count_frip(tags,self._peaks("nuclear\t0\t100\n")),(2,2,1.0))

    def _peaks(self,content):
        path=self.d/("peaks_"+str(len(list(self.d.glob("peaks_*"))))+".bed")
        path.write_text(content)
        return path

    def test_empty_nuclear_reference_fails(self):
        with self.assertRaises(ValueError):nuclear_references(["NC_007936.1","*"],"NC_007936.1")

    def test_nuclear_span_is_computed_and_excludes_mitochondria(self):
        fai=self.d/"reference.fai"
        fai.write_text("chr1\t100\t0\t0\t0\nNC_007936.1\t20\t0\t0\t0\nchr2\t50\t0\t0\t0\n")
        self.assertEqual(calculate_reference_span(fai,"NC_007936.1"),([("chr1",100),("chr2",50)],150))

    def test_tn5_plus_minus_and_boundaries(self):
        self.assertEqual(insertion_position(100,130,False),104)
        self.assertEqual(insertion_position(100,130,True),125)
        self.assertEqual(insertion_position(0,3,True),-2)

    def test_tss_boundaries_and_malformed(self):
        bed=self.d/"tss.bed";bed.write_text("chr1\t0\t1\ta\t0\t+\nchr1\t20\t21\tb\t0\t-\n")
        result=load_tss_bed(bed,{"chr1":40},10)
        self.assertEqual(result[2:],(2,1,1,0))
        bed.write_text("chr1\t10\t12\ta\t0\t+\n")
        with self.assertRaises(ValueError):load_tss_bed(bed,{"chr1":40},10)
        bed.write_text("chr1\t10\t11\ta\t0\t?\n")
        with self.assertRaises(ValueError):load_tss_bed(bed,{"chr1":40},10)

    def test_tss_profile_validation_rejects_malformed_data(self):
        profile=self.d/"profile.tsv";summary=self.d/"summary.tsv"
        profile.write_text("relative_start\trelative_end\tinsertion_count\tnormalized_signal\n-10\t0\t1\t1.0\n0\t10\t2\t2.0\n")
        summary.write_text("metric\tvalue\nusable_tss\t1\nwindow_bp\t10\nbin_size_bp\t10\nbackground_bp_per_side\t10\nbackground_mean\t1\ntss_enrichment_score\t2\n")
        validate_tss(profile,summary,10,10,10)
        profile.write_text(profile.read_text().replace("1.0","nan"))
        with self.assertRaises(ValueError):validate_tss(profile,summary,10,10,10)

    def test_tss_rejects_alignment_flags_and_low_mapq(self):
        class TSSRead:
            def __init__(self,flag="",mapq=30):
                self.is_unmapped=flag=="unmapped";self.is_secondary=flag=="secondary"
                self.is_supplementary=flag=="supplementary";self.is_duplicate=flag=="duplicate"
                self.is_qcfail=flag=="qcfail";self.mapping_quality=mapq
        self.assertTrue(usable_read(TSSRead(),30))
        for flag in ("unmapped","secondary","supplementary","duplicate","qcfail"):
            self.assertFalse(usable_read(TSSRead(flag),30))
        self.assertFalse(usable_read(TSSRead(mapq=29),30))

    def test_frip_deterministic(self):
        tags=self.d/"tags.bed";peaks=self.d/"peaks.bed"
        tags.write_text("chr1\t1\t5\nchr1\t8\t10\nchr1\t20\t21\nchr2\t1\t2\n")
        peaks.write_text("chr1\t4\t9\nchr2\t5\t6\n")
        self.assertEqual(count_frip(tags,peaks),(4,2,0.5))

    def test_bigwig_boundary_clipping_and_unknown_rejection(self):
        sizes=self.d/"sizes";sizes.write_text("chr1\t100\n")
        bg=self.d/"x.bdg";out=self.d/"out.bdg";metrics=self.d/"m.json"
        bg.write_text("chr1\t-5\t10\t1\nchr1\t90\t110\t2\nchr1\t120\t130\t3\n")
        value=clip(bg,sizes,out,metrics)
        self.assertEqual((value["input_intervals"],value["kept_intervals"],value["clipped_intervals"],value["dropped_intervals"],value["clipped_bp_total"]),(3,2,3,1,45))
        bg.write_text("unknown\t0\t1\t1\n")
        with self.assertRaises(ValueError):clip(bg,sizes,out,metrics)

    def test_bigwig_valid_signal_subset_and_unused_reference_allowed(self):
        validate_header({"chr1":100},{"chr1":100,"chr2":200},{"chr1"})

    def test_bigwig_file_contract_valid_subset(self):
        import pyBigWig
        sizes=self.d/"header.sizes";sizes.write_text("chr1\t100\nchr2\t200\n")
        signal=self.d/"header.bedGraph";signal.write_text("chr1\t0\t5\t1\n")
        bigwig=self.d/"header.bw"
        with pyBigWig.open(str(bigwig),"w") as handle:
            handle.addHeader([("chr1",100)])
            handle.addEntries(["chr1"],[0],ends=[5],values=[1.0])
        validate_bigwig(bigwig,sizes,signal)

    def test_bigwig_missing_signal_chromosome_fails(self):
        with self.assertRaisesRegex(ValueError,"missing signal"):
            validate_header({"chr1":100},{"chr1":100,"chr2":200},{"chr1","chr2"})

    def test_bigwig_extra_chromosome_fails(self):
        with self.assertRaisesRegex(ValueError,"outside the reference"):
            validate_header({"chr1":100,"chrX":50},{"chr1":100},{"chr1"})

    def test_bigwig_reference_chromosome_without_signal_in_header_fails(self):
        with self.assertRaisesRegex(ValueError,"without retained signal"):
            validate_header({"chr1":100,"chr2":200},{"chr1":100,"chr2":200},{"chr1"})

    def test_bigwig_chromosome_length_mismatch_fails(self):
        with self.assertRaisesRegex(ValueError,"length mismatch"):
            validate_header({"chr1":99},{"chr1":100},{"chr1"})

    def test_bigwig_signal_chromosomes_are_read_from_clipped_bedgraph(self):
        bedgraph=self.d/"signal.bedGraph"
        bedgraph.write_text("chr1\t0\t5\t1\nchr3\t5\t8\t2\n")
        self.assertEqual(signal_chromosomes(bedgraph),{"chr1","chr3"})

    def test_streaming_sha256_helper(self):
        path=self.d/"hash.bin";content=(b"abcdef0123456789"*1024)+b"tail"
        path.write_bytes(content)
        self.assertEqual(sha256_file(path,chunk_size=7),hashlib.sha256(content).hexdigest())
        with self.assertRaises(ValueError):sha256_file(path,chunk_size=0)

    def test_replicate_group_and_idr_eligibility(self):
        groups=load_groups(ROOT/"config/atacseq_replicate_groups.tsv",self.manifest)
        self.assertEqual(groups["PRJNA667472_CHO-K1_bulk_ATAC"],["SRR12774931","SRR12774932"])
        bad=self.d/"groups.tsv";bad.write_text("replicate_group\trun_accession\tcompatibility\treplicate_type\nG\tSRR29929613\tcompatible\tunspecified\nG\tSRR12774931\tcompatible\tunspecified\n")
        with self.assertRaises(ValueError):load_groups(bad,self.manifest)

    def test_idr_preserves_both_thresholds(self):
        source=self.d/"idr.tsv"
        def row(score):return "\t".join(["chr1","0","1","p","1",".","1","1","1","1","1",str(score)])+"\n"
        source.write_text(row(0.9)+row(1.0)+row(1.30103))
        out10=self.d/"10.tsv";out05=self.d/"05.tsv";metrics=self.d/"idr.json"
        filter_idr(source,out10,out05,metrics)
        self.assertEqual(len(out10.read_text().splitlines()),2)
        self.assertEqual(len(out05.read_text().splitlines()),1)

    def test_real_production_snakemake_dag_direct_downstream_gate(self):
        snakemake=shutil.which("snakemake")
        self.assertIsNotNone(snakemake,"snakemake is required for the ATAC DAG regression")
        provenance=self.d/"atac.provenance.json";provenance.write_text("{}\n")
        config=self.d/"dryrun.yaml"
        config.write_text(
          "fastq:\n  manifest: "+str(self.source_manifest)+"\n"
          "atacseq:\n  production_manifest: "+str(self.manifest)+"\n"
          "  production_manifest_provenance: "+str(provenance)+"\n")
        environment=os.environ.copy();environment["XDG_CACHE_HOME"]=str(self.d/"cache")
        base=[snakemake,"--snakefile",str(ROOT/"Snakefile"),"--dry-run","--cores","1"]
        config_args=["--configfile",str(config),"--default-resources","tmpdir=/tmp"]
        for run in ("SRR12774931","SRR12774932"):
            target=f"results/atacseq/provenance/{run}/provenance.json"
            result=subprocess.run(base+[target]+config_args,text=True,capture_output=True,env=environment,cwd=ROOT)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertIn("atacseq_run_provenance",result.stdout+result.stderr)
        for run in ("SRR12774934","SRR29929613"):
            target=f"results/atacseq/provenance/{run}/provenance.json"
            result=subprocess.run(base+[target]+config_args,text=True,capture_output=True,env=environment,cwd=ROOT)
            self.assertNotEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertIn(f"ATAC run {run} is not production eligible",result.stdout+result.stderr)


if __name__=="__main__":unittest.main()
