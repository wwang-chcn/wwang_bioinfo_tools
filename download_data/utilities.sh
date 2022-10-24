#! /bin/bash

function ena_meta_download {
    project_accession=${1}
    wget -O ${project_accession}.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${project_accession}&result=read_run&fields=sample_accession,experiment_accession,run_accession,scientific_name,library_layout,fastq_aspera,fastq_ftp&format=tsv&download=true"
}