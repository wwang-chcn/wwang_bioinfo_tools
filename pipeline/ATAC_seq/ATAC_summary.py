#!/usr/bin/env python

"""
wrote by Wen Wang (wwang.tongji@gmail.com)
"""

# Re-write: Dex-12-2020
# Re-write: July-19-2019

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
            _, _, name, processer, genomeVersion, reads1_files, reads2_files, *_ = line.strip().split()
            yield name, genomeVersion, reads1_files, reads2_files


def source_files_check(name, genomeVersion, reads1_files, reads2_files):
    '''source files for pipeline running check'''
    flag = True

    # mapping index
    # XX.[1-4].bt2 XX.rev.[12].bt2
    index_files_suffix = ['.1.bt2','.2.bt2','.3.bt2','.4.bt2','.rev.1.bt2','.rev.2.bt2']
    for suffix in index_files_suffix:
        if not os.path.isfile(os.path.expanduser(f'~/source/bySpecies/{genomeVersion}/{genomeVersion}{suffix}')):
            os.sys.stdout(f'Error! No bowtie2 index file: {genomeVersion}{suffix} in {os.path.expanduser(f"~/source/bySpecies/{genomeVersion}")}\n')
            flag = False

    # read_file_check
    reads1_files = reads1_files.split(',')
    reads2_files = reads2_files.split(',')
    if not len(reads1_files) == len(reads2_files):
        os.sys.stdout(f'Error! Number of input reads file1 and reads file2 do not match!\n')
        flag = False
    for i in range(len(reads1_files))
        if reads1_files[i] == reads2_files[i]:
            os.sys.stdout(f'Error! The {i}-th input reads file1 and reads file2 {reads1_files[i]} were same!\n')
            flag = False
        if not os.path.isfile(f'./0_raw_data/{reads1_files[i]}'):
            os.sys.stdout(f'Error! Input reads file {reads1_files[i]} do not exist!\n')
            flag = False
        if not os.path.isfile(f'./0_raw_data/{reads2_files[i]}'):
            os.sys.stdout(f'Error! Input reads file {reads2_files[i]} do not exist!\n')
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
    if len(sys.argv) < 3:
        sys.stdout.write('No enought arguments!\n')
        sys.stdout.write(USAGE+'\n')
        sys.exit(1)
    
    outputPrefix = sys.argv[1]
    
    with open('ATAC_summary_{}_mapping_info.csv'.format(outputPrefix),'w') as mappingFhd, \
         open('ATAC_summary_{}_fragments_length.csv'.format(outputPrefix), 'w') as flFhd:
        mapping_csv = csv.writer(mappingFhd)
        fl_rows = []
        # header
        mapping_csv.writerow(['Sample','Label','Raw Reads Pair (adapter filtered)','Mapped Reads Pair (mapped cordantly, q30 filtered)','Mapping Efficiency','chrM Fragments','chrM Percentage','Unique Nuclear Fragments','Duplicate Level','Effective Fragments'])
        
        for name, genomeVersion, reads1_files, reads2_files in load_sample('runned.sh'):
            sys.stdout.write(line) # for debugging
            if not source_files_check(name, genomeVersion, reads1_files, reads2_files):
                continue
            sample, label = name, name # running name, presenting label
            output_line = [sample,label] # sample \t label
            
            # mapping & chrosome M information
            total, mapped, chrM = 0, 0, 0
            with open('./1_mapping/{}_mapping.log'.format(sample), 'r') as fhd:
                line = next(fhd)
                total = int(line.strip().split()[0])
            with open('./2_signal/{}_chromosome_distribution.txt'.format(sample), 'r') as cd:
                cd.readline()
                for line in cd:
                    line = line.strip().split()
                    if line[0] == 'chrM':
                        chrM += int(line[1])
                    mapped += int(line[1])
            output_line.extend([total, mapped, 1.0 * mapped / total, chrM, 1.0 * chrM / mapped]) # Raw reads pair \t mapped reads pair \t mapping efficiency \t chrM \t chrM percentage
            
            # unique nuclear fragments & duplicate level
            raw_fragment = int(subprocess.check_output(f'''bigBedInfo 2_signal/{sample}_raw_fragments.bb | grep itemCount | cut -d " " -f 2 | sed -s 's/,//g' ''', shell=True).decode().strip())
            fragment = int(subprocess.check_output(f'''bigBedInfo 2_signal/{sample}_fragments.bb | grep itemCount | cut -d " " -f 2 | sed -s 's/,//g' ''', shell=True).decode().strip())
            output_line.extend([fragment, 1 - 1.0 * fragment / raw_fragment]) # unique nuclear fragments \t duplicate level

            OCR_fragment = int(subprocess.check_output(f'''bigBedInfo 2_signal/OCR/{sample}_uniq_OCR_SE_reads.bb | grep itemCount | cut -d " " -f 2 | sed -s 's/,//g' ''', shell=True).decode().strip())
            output_line.append(OCR_fragment) # Effective Fragments
            
            mapping_csv.writerow(output_line)
            
            with open('./2_signal/{}_fragments_length.txt'.format(sample), 'r') as fl:
                next(fl)
                fL = [0 for i in range(301)]
                for line in fl:
                    line = line.strip().split()
                    try:
                        fL[int(line[0])] = int(line[1])
                    except IndexError:
                        pass
                fL.insert(0,sample)
            fl_rows.append(fL)
        headers = ['Sample']
        headers.extend(range(301))
        fl_csv = csv.writer(flFhd)
        fl_csv.writerow(headers)
        fl_csv.writerows(fl_rows)


# ------------------------------------
# Program running
# ------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you ^.^!")
        sys.exit(0)




