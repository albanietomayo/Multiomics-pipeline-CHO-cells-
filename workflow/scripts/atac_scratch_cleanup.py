import shutil
import sys
from pathlib import Path
RUNS = {"SRR12774931", "SRR12774932"}

def clean(run, work, scratch):
    root = Path(work).resolve()
    tmp = Path(scratch).resolve()
    if run not in RUNS or root.parent != tmp or Path.cwd().resolve() != root or not root.name.startswith('atac_bulk_'):
        raise ValueError('Scratch cleanup scope mismatch')
    for relative in ('data/raw/fastq','results/preprocessing/fastp','results/atacseq/bowtie2','results/atacseq/filtered_bam'):
        path = root / relative / run
        if path.is_symlink() or root not in path.resolve().parents:
            raise ValueError('Unsafe scratch path')
        if path.exists():
            shutil.rmtree(path)

if __name__ == '__main__': clean(*sys.argv[1:])
