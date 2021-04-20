# Dependencies

## Softwares

- Python (≥3.6)
- Cutadapt (≥3.3)
- FastQC (≥0.11.4)
- TrimGalore (≥0.6.4)
- Bowtie2 (≥2.4.2)
- htslib (≥1.6)
- SAMtools (≥1.6)
- BEDtools (≥2.27.1)
- GNU Awk (≥4.1.3)

## Input & source files

- Raw sequencing data (paired-end) in ```0_raw_data```
- Bowtie2 mapping index in ```~/source/bySpecies/${genomeVersion}/``` with bt2_index_base name of ```${genomeVersion}``` (whole genome) and ```${genomeVersion}_main``` (genome without small scaffold)
- Chrom size files (```${genomeVersion}.chorm.sizes ```and ```${genomeVersion}_main.chorm.sizes```)

# Output files

Files in ```2_signal/OCR``` contain the chromatin accessibility information.