# Dependencies

- [Moabs](https://github.com/sunnyisgalaxy/moabs) (bsmap, mcall)

- [Trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [samtools](http://www.htslib.org/)

# Usage

```bash
DNAme_PE.sh <name> <processer> <sourceDir> <genomeVersion> <reads1,+> <reads2,+>
```

## Example

```bash
bash DNAme_PE.sh test 32 sourceDir mm10 test_r1_1.fastq.gz,test_r2_1.fastq.gz,test_r3_1.fastq.gz,test_r4_1.fastq.gz test_r1_2.fastq.gz,test_r2_2.fastq.gz,test_r3_2.fastq.gz,test_r4_2.fastq.gz
```

# Files required in source dir

genomeVersion_main.fa (scaffold and small chromsome removed. You can generate this file using [getMainFa.py](https://github.com/wwang-chcn/wwang_bioinfo_tools/blob/master/genome_source_build/getMainFa.py))

# FAQ

