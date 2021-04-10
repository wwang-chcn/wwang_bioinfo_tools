# Dependencies

## Softwares

- Python2 (2.7) & Python3 ($\geq$ 3.6)
- Cutadapt ($\geq$ 3.3)
- FastQC ($\geq$ 0.11.4)
- TrimGalore ($\geq$ 0.6.4)
- Bowtie2 ($\geq$ 2.4.2)
- htslib ($\geq$ 1.6)
- SAMtools ($\geq$ 1.6)
- BEDtools ($\geq$ 2.27.1)
- GNU Awk ($\geq$ 4.1.3)
- MACS2 ($\geq$ 2.1.1.20160309)

## Input & source files

- Raw sequencing data in ```0_raw_data```
- Bowtie2 mapping index in ```~/source/bySpecies/${genomeVersion}/``` with bt2_index_base name of ```${genomeVersion}``` (whole genome) and ```${genomeVersion}_main``` (genome without small scaffold)
- Chrom size files (```${genomeVersion}.chorm.sizes ```and ```${genomeVersion}_main.chorm.sizes```)

