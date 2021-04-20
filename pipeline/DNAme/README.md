# Dependencies

- [Moabs](https://github.com/sunnyisgalaxy/moabs) (bsmap, mcall, ≥1.3.9.6)
- Cutadapt (≥3.3)
- FastQC (≥0.11.4)
- TrimGalore (≥0.6.4)
- htslib (≥1.6)
- SAMtools (≥1.6)

# Usage

```bash
DNAme_PE.sh <name> <processer> <sourceDir> <genomeVersion> <reads1,+> <reads2,+>
```

## Example

```bash
bash DNAme_PE.sh test 32 sourceDir mm10 test_r1_1.fastq.gz,test_r2_1.fastq.gz,test_r3_1.fastq.gz,test_r4_1.fastq.gz test_r1_2.fastq.gz,test_r2_2.fastq.gz,test_r3_2.fastq.gz,test_r4_2.fastq.gz
```

## Input & source files

- Raw sequencing data (paired-end) in ```0_raw_data```
- Genome sequence file ```${genomeVersion}_main.fa``` (fastq format) in ```~/source/bySpecies/${genomeVersion}/```

