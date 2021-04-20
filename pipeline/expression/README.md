# Dependencies

## Softwares

- Python2 (2.7) & Python3 (≥3.6)
- Cutadapt (≥3.3)
- FastQC (≥0.11.4)
- TrimGalore (≥0.6.4)
- HISAT2 (≥2.1.0)
- StringTie (≥1.3.3b)
- GFOLD (≥1.1.4)
- htslib (≥1.6)
- SAMtools (≥1.6)
- BEDtools (≥2.27.1)
- RSeQC (≥2.6.4)

## Input & source files

- Raw sequencing data in ```0_raw_data```
- HISAT2 mapping index in ```~/source/bySpecies/${genomeVersion}/``` with ht2_index_base name of ```${genomeVersion}_main``` (genome without small scaffold)
- refGene annotation file ```${genomeVersion}.refGene.gtf``` in ```~/source/bySpecies/${genomeVersion}/```
- ensGene annotation file ```${genomeVersion}.ensGene.gtf``` in ```~/source/bySpecies/${genomeVersion}/``` (optional)
- Chrom size files (```${genomeVersion}.chorm.sizes ```and ```${genomeVersion}_main.chorm.sizes```)

