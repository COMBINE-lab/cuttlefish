# E. coli colored benchmark

This fixture pins 150 complete, non-atypical RefSeq assemblies returned by
NCBI Datasets for `Escherichia coli`. Each assembly is a separate graph color.

Prepare the genomic FASTA files with:

```bash
scripts/prepare_ecoli_color_benchmark.sh
```

The command writes `genomes.list`, with one FASTA path per source. Downloaded
sequence data is intentionally excluded from Git. `assemblies.tsv` records the
assembly metadata returned with the pinned accessions.
