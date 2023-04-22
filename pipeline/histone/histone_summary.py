#!/usr/bin/env python3

# ------------------------------------
# Load Modules
# ------------------------------------
import os
import sys
import csv
import subprocess


# ------------------------------------
# Sub Functions
# ------------------------------------
def load_sample(running_file):
    with open(running_file) as fhd:
        for line in fhd:
            if not line.startswith('bash'):
                continue
            _, _, name, processer, genomeVersion, reads1_files, reads2_files, *_ = line.strip(
            ).split()
            yield name, genomeVersion, reads1_files, reads2_files


def source_files_check(name, genomeVersion, reads1_files, reads2_files):
    '''source files for pipeline running check'''
    flag = True

    # mapping index
    # XX.[1-4].bt2 XX.rev.[12].bt2
    index_files_suffix = [
        '.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2'
    ]
    for suffix in index_files_suffix:
        if not os.path.isfile(
                os.path.expanduser(
                    f'~/source/bySpecies/{genomeVersion}/{genomeVersion}{suffix}'
                )):
            os.sys.stdout.write(
                f'Error! No bowtie2 index file: {genomeVersion}{suffix} in {os.path.expanduser(f"~/source/bySpecies/{genomeVersion}")}\n'
            )
            flag = False

    # read_file_check
    reads1_files = reads1_files.split(',')
    reads2_files = reads2_files.split(',')
    if len(reads1_files) != len(reads2_files):
        os.sys.stdout.write(
            'Error! Number of input reads file1 and reads file2 do not match!\n'
        )
        flag = False
    for i in range(len(reads1_files)):
        if reads1_files[i] == reads2_files[i]:
            os.sys.stdout.write(
                f'Warning! The {i}-th input reads file1 and reads file2 {reads1_files[i]} were same!\n'
            )
            flag = False
        if not os.path.isfile(f'./0_raw_data/{reads1_files[i]}'):
            os.sys.stdout.write(
                f'Warning! Input reads file {reads1_files[i]} do not exist!\n')
            flag = False
        if not os.path.isfile(f'./0_raw_data/{reads2_files[i]}'):
            os.sys.stdout.write(
                f'Warning! Input reads file {reads2_files[i]} do not exist!\n')
            flag = False

    return flag


def check_mapping():
    pass


# ------------------------------------
# Main Function
# ------------------------------------


def main():

    script_name = os.path.basename(sys.argv[0])
    USAGE = f'{script_name} <outputPrefix>'
    if len(sys.argv) < 2:
        sys.stdout.write('No enought arguments!\n')
        sys.stdout.write(USAGE + '\n')
        sys.exit(1)

    outputPrefix = sys.argv[1]

    with open(f'Histone_summary_{outputPrefix}_mapping_info.csv', 'w') as mappingFhd:
        mapping_csv = csv.writer(mappingFhd)
        fl_rows = []
        # header
        mapping_csv.writerow([
            'Sample', 'Label', 'Raw Reads Pair (adapter filtered)',
            'Mapped Reads (Pair) (mapped cordantly, q30 filtered)',
            'Mapping Efficiency', 'Unique fragments', 'Duplicate level', 'Peak number'
        ])

        for name, genomeVersion, reads1_files, reads2_files in load_sample(
                'runned.sh'):
            sys.stdout.write(
                f'name: {name}\n\tgenome version: {genomeVersion}\n\treads1 files: {reads1_files}\n\treads2 files: {reads2_files}\n'
            )  # for debugging
            if not source_files_check(name, genomeVersion, reads1_files,
                                      reads2_files):
                pass
            sample, label = name, name  # running name, presenting label
            output_line = [sample, label]  # sample \t label

            # mapping 
            with open(f'./1_mapping/{sample}_mapping.log', 'r') as fhd:
                line = next(fhd)
                total = int(line.strip().split()[0])
            if os.path.isfile(f'./2_signal/{name}_raw_fragments.bed'):
                mapped = int(subprocess.check_output(f'cat ./2_signal/{name}_raw_fragments.bed | wc -l', shell=True).decode().strip())
            elif os.path.isfile(f'./2_signal/{name}_raw_fragments.bb'):
                mapped = subprocess.check_output(f'''bigBedInfo ./2_signal/{name}_raw_fragments.bb | grep itemCount | cut -d ' ' -f 2''', shell=True).decode().strip()
                mapped = int(mapped.replace(',',''))
            elif os.path.isfile(f'./2_signal/{name}_raw_reads.bed'):
                mapped = int(subprocess.check_output(f'cat ./2_signal/{name}_raw_reads.bed | wc -l', shell=True).decode().strip())
            elif os.path.isfile(f'./2_signal/{name}_raw_reads.bb'):
                mapped = subprocess.check_output(f'''bigBedInfo ./2_signal/{name}_raw_reads.bb | grep itemCount | cut -d ' ' -f 2''', shell=True).decode().strip()
                mapped = int(mapped.replace(',',''))
            else:
                sys.stderr.write(f'No raw fragments/reads file found for sample: {name}\n')
                continue
            
            # unique fragments/reads & duplicate level
            if os.path.isfile(f'./2_signal/{name}_fragments.bed'):
                unique = int(subprocess.check_output(f'cat ./2_signal/{name}_fragments.bed | wc -l', shell=True).decode().strip())
            elif os.path.isfile(f'./2_signal/{name}_fragments.bb'):
                unique = subprocess.check_output(f'''bigBedInfo ./2_signal/{name}_fragments.bb | grep itemCount | cut -d ' ' -f 2''', shell=True).decode().strip()
                unique = int(unique.replace(',',''))
            elif os.path.isfile(f'./2_signal/{name}_reads.bed'):
                unique = int(subprocess.check_output(f'cat ./2_signal/{name}_reads.bed | wc -l', shell=True).decode().strip())
            elif os.path.isfile(f'./2_signal/{name}_reads.bb'):
                unique = subprocess.check_output(f'''bigBedInfo ./2_signal/{name}_reads.bb | grep itemCount | cut -d ' ' -f 2''', shell=True).decode().strip()
                unique = int(unique.replace(',',''))
            else:
                sys.stderr.write(f'No fragments/reads file found for sample: {name}\n')
                continue

            # peak number
            if os.path.isfile(f'./3_peak/{name}_peaks.narrowPeak'):
                peak_num = int(subprocess.check_output(f'cat ./3_peak/{name}_peaks.narrowPeak | wc -l', shell=True).decode().strip())
            elif os.path.isfile(f'./3_peak/{name}_peaks.broadPeak'):
                peak_num = int(subprocess.check_output(f'cat ./3_peak/{name}_peaks.broadPeak | wc -l', shell=True).decode().strip())
            else:
                sys.stderr.write(f'No peaks file found for sample: {name}\n')
                continue
            output_line.extend(
                [
                    total, mapped, mapped / total, unique, 1 - unique / mapped, peak_num
                ]
            )  # Raw reads pair \t mapped reads pair \t mapping efficiency 



            mapping_csv.writerow(output_line)
        headers = ['Sample']
        headers.extend(range(301))


# ------------------------------------
# Program running
# ------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stdout.write("User interrupts me! ;-) See you ^.^!")
        sys.exit(0)
