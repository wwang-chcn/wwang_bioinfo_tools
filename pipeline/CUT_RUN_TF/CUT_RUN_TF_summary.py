#!/usr/bin/env python3

# ------------------------------------
# Load Modules
# ------------------------------------
import os
import sys
import csv
import subprocess

from typing import Optional


# ------------------------------------
# Sub Functions
# ------------------------------------
class MappingStatus(object):
    def __init__(self) -> None:
        self.total: Optional[int] = None
        self.mapped: Optional[int] = None
        self.chrM: Optional[int] = None
        self.status = False


class FrgmentsStatus(object):
    def __init__(self) -> None:
        self.raw: Optional[int] = None
        self.uniq: Optional[int] = None
        self.OCR: Optional[int] = None
        self.status = True


class PeakStatus(object):
    def __init__(self) -> None:
        self.peak: Optional[int] = None
        self.FFiP: Optional[float] = None
        self.status = True

class CRTFsample(object):
    def __init__(self, name: str, genomeVersion: str, reads1_files: str, reads2_files: str) -> None:
        self.name = name
        self.genomeVersion = genomeVersion
        self.reads1_files = reads1_files.split(',')
        self.reads2_files = reads2_files.split(',')
        self.mapping_status = MappingStatus()
        self.fragments_status = FrgmentsStatus()
        self.peak_status = PeakStatus()
        self.OCR_peak_status = PeakStatus()
    
    def file_check(self, file: str, file_type: str) -> bool:
        """check file existence related to this sample"""
        if os.path.isfile(file):
            return True
        else:
            sys.stdout.write(f'Error! The {file_type} file ({file}) of sample: {self.name} do not exist! Skipped.\n')
            return False

    def source_files_check(self) -> bool:
        """source files for pipeline running check"""
        flag = True
        # mapping index check
        # XX.[1-4].bt2 XX.rev.[12].bt2
        index_files_suffix = [
            '.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2'
        ]
        for suffix in index_files_suffix:
            if not os.path.isfile(
                    os.path.expanduser(
                        f'~/source/bySpecies/{self.genomeVersion}/{self.genomeVersion}{suffix}'
                    )):
                sys.stdout.write(
                    f'Error! No bowtie2 index file: {self.genomeVersion}{suffix} in {os.path.expanduser(f"~/source/bySpecies/{self.genomeVersion}")}\n'
                )
                flag = False
        # reads file check
        if len(self.reads1_files) != len(self.reads2_files):
            sys.stdout.write(
                'Error! Number of input reads file1 and reads file2 do not match!\n'
            )
            flag = False
        for i in range(len(self.reads1_files)):
            if self.reads1_files[i] == self.reads2_files[i]:
                sys.stdout.write(
                    f'Warning! The {i}-th input reads file1 and reads file2 {self.reads1_files[i]} were same!\n'
                )
                flag = False
            if not os.path.isfile(f'./0_raw_data/{self.reads1_files[i]}'):
                sys.stdout.write(
                    f'Warning! Input reads file {self.reads1_files[i]} do not exist!\n')
                flag = False
            if not os.path.isfile(f'./0_raw_data/{self.reads2_files[i]}'):
                sys.stdout.write(
                    f'Warning! Input reads file {self.reads2_files[i]} do not exist!\n')
                flag = False
        return flag
    
    def get_input_reads_info(self) -> None:
        """get the number of input reads pairs from mapping log file"""
        mapping_log_file = f'1_mapping/Mapping_{self.name}.log'
        # mapping log file check
        if not self.file_check(mapping_log_file, 'mapping log'):
            return
        # total reads pairs
        with open(mapping_log_file, 'r') as fhd:
                line = next(fhd)
                self.mapping_status.total = int(line.strip().split()[0])
    
    def get_mapped_info(self) -> None:
        """get total & chrM mapped fragments num from chromsome distribution file"""
        chromsome_distribution_file = f'4_basic_QC/{self.name}_raw_chromosome_distribution.txt'
        # chromsome distribution file check
        if not self.file_check(chromsome_distribution_file, 'chromsome distribution'):
            return
        # get total & chrM mapped fragments num
        mapped, chrM = 0, 0
        with open(chromsome_distribution_file, 'r') as fhd:
                fhd.readline()
                for line in fhd:
                    line = line.strip().split()
                    if line[0] == 'chrM':
                        chrM += int(line[1])
                    mapped += int(line[1])
        self.mapping_status.mapped = mapped
        self.mapping_status.chrM = chrM
        self.mapping_status.status = True
    
    def get_num_of_bed_file(self, bed_file: str) -> int:  #TODO: multiple type hint
        bed_file_num = 0
        bigbed_file = f'{bed_file[:-2]}b'
        if self.file_check(f'{bigbed_file}', 'BigBed file'):
            bed_file_num = int(
                subprocess.check_output(
                    f'''bigBedInfo {bigbed_file} | grep itemCount | cut -d " " -f 2 | sed -s 's/,//g' ''',
                    shell=True).decode().strip())
            return bed_file_num
        elif self.file_check(f'{bed_file}', 'Bed file'):
            bed_file_num = int(
                subprocess.check_output(
                    f'''cat {bed_file} | wc -l''',
                    shell=True).decode().strip())
            return bed_file_num
        else:
            sys.stdout.write(f'Warning! Raw fragments file ({bigbed_file} or {bed_file}) of sample: {self.name} do not exist, skip this sample.\n')
            raise FileNotFoundError  #TODO: return None when using nultiple type hint
    
    def get_fragments_num(self) -> None:
        """get raw, uniq, OCR fragments num"""
        # raw fragments
        try:
            self.fragments_status.raw = self.get_num_of_bed_file(f'2_signal/{self.name}_raw_fragments.bed')
        except FileNotFoundError:
            self.fragments_status.status = False
        # uniq fragments
        try:
            self.fragments_status.uniq = self.get_num_of_bed_file(f'2_signal/{self.name}_fragments.bed')
        except FileNotFoundError:
            self.fragments_status.status = False
        # raw fragments
        try:
            self.fragments_status.OCR = self.get_num_of_bed_file(f'2_signal/{self.name}_OCR_fragments.bed')
        except FileNotFoundError:
            self.fragments_status.status = False
    
    def cutsites_files_check(self) -> None:
        """cut sites files check"""
        if not os.path.isfile(f'2_signal/{self.name}_cutsites.bed') and not os.path.isfile(f'2_signal/{self.name}_cutsites.bb'):
            sys.stdout.write(f'Warning! Cut sites files: ({self.name}_cutsites.bed and {self.name}_cutsites.bb) for sample: {self.name} were not exist!\n')
        if not os.path.isfile(f'2_signal/{self.name}_cutsites_plus.bed') and not os.path.isfile(f'2_signal/{self.name}_cutsites_plus.bb'):
            sys.stdout.write(f'Warning! Cut sites plus strand files: ({self.name}_cutsites_plus.bed and {self.name}_cutsites_plus.bb) for sample: {self.name} were not exist!\n')
        if not os.path.isfile(f'2_signal/{self.name}_cutsites_minus.bed') and not os.path.isfile(f'2_signal/{self.name}_cutsites_minus.bb'):
            sys.stdout.write(f'Warning! Cut sites minus strand files: ({self.name}_cutsites_minus.bed and {self.name}_cutsites_minus.bb) for sample: {self.name} were not exist!\n')

    def get_number_of_peak_file(self, peak_file: str) -> int:
        if self.file_check(peak_file, 'Peak file'):
            peak_num = int(subprocess.check_output(
                    f'''cat {peak_file} | wc -l''',
                    shell=True).decode().strip())
            return peak_num

    def get_peak_info(self) -> None:
        """peak info (all fragments by default)"""

        peak_file = f'3_peak/{self.name}_peaks.narrowPeak'

        try:
            self.peak_status.peak = self.get_number_of_peak_file(peak_file)
        except FileNotFoundError:
            self.peak_status.status = False

        fragments_file = f'2_signal/{self.name}_fragments.bed'
        fragments_bb_file = fragments_file[:-2] + 'b'
        if not os.path.isfile(fragments_file):
            if os.path.isfile(fragments_bb_file):
                subprocess.Popen(f'bigBedToBed {fragments_bb_file} {fragments_file}', shell=True)
            else:
                raise FileNotFoundError(f'Fragments file ({fragments_file}) not found')
        FiP = int(subprocess.check_output(
                    f'''intersectBed -u -a {fragments_file} -b {peak_file} | wc -l''',
                    shell=True).decode().strip())
        self.peak_status.FFiP = 1.0 * FiP / self.fragments_status.uniq

    def get_OCR_peak_info(self) -> None:
        """peak info (OCR fragments)"""

        peak_file = f'3_peak/{self.name}_OCR_peaks.narrowPeak'

        try:
            self.OCR_peak_status.peak = self.get_number_of_peak_file(peak_file)
        except FileNotFoundError:
            self.OCR_peak_status.status = False

        fragments_file = f'2_signal/{self.name}_OCR_fragments.bed'
        fragments_bb_file = fragments_file[:-2] + 'b'
        if not os.path.isfile(fragments_file):
            if os.path.isfile(fragments_bb_file):
                subprocess.Popen(f'bigBedToBed {fragments_bb_file} {fragments_file}', shell=True)
            else:
                raise FileNotFoundError(f'Fragments file ({fragments_file}) not found')
        FiP = int(subprocess.check_output(
                    f'''intersectBed -u -a {fragments_file} -b {peak_file} | wc -l''',
                    shell=True).decode().strip())
        self.OCR_peak_status.FFiP = 1.0 * FiP / self.fragments_status.uniq
    
    def output_line(self) -> list:
        """output in list format"""
        out = []
        if not all([self.mapping_status.status, self.fragments_status.status]):
            sys.stdout.write(f'Warning! Infomation about sample: {self.name} has not enought. Skipped.\n')
            return out
        out.extend([self.name, self.genomeVersion])
        out.extend([self.mapping_status.total, self.mapping_status.mapped, self.mapping_status.mapped / self.mapping_status.total, self.mapping_status.chrM])
        out.extend([self.fragments_status.raw, self.fragments_status.uniq, 1 - self.fragments_status.uniq / self.fragments_status.raw, self.fragments_status.OCR])
        out.extend([self.peak_status.peak, self.peak_status.FFiP, self.OCR_peak_status.peak, self.OCR_peak_status.FFiP])
        return out


# ------------------------------------
# Sub Functions
# ------------------------------------
def load_sample(running_file):
    with open(running_file) as fhd:
        for line in fhd:
            if not line.startswith('bash'):
                continue
            _, _, name, _, genomeVersion, reads1_files, reads2_files, *_ = line.strip(
            ).split()
            yield name, genomeVersion, reads1_files, reads2_files


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

    with open(f'CRTF_summary_{outputPrefix}_mapping_info.csv', 'w') as mappingFhd, \
         open(f'CRTF_summary_{outputPrefix}_fragments_length.csv', 'w') as flFhd:
        mapping_csv = csv.writer(mappingFhd)
        fl_rows = []
        # header
        mapping_csv.writerow([
            'Sample', 'Genome assembly', 'Raw reads pair (adapter removed)',
            'Mapped reads pair (mapped cordantly, q30 filtered)',
            'Mapping efficiency', 'chrM fragments', 'Raw fragments',
            'Unique nuclear fragments', 'Duplicate level', 'Effective fragments',
            'Peak number (all fragments)', 'FFiP (all fragments)', 'Peak number (OCR fragments)', 'FFiP (OCR fragments)'
        ])

        for name, genomeVersion, reads1_files, reads2_files in load_sample(
                'runned.sh'):
            # sys.stdout.write(
            #     f'name: {name}\n\tgenome version: {genomeVersion}\n\treads1 files: {reads1_files}\n\treads2 files: {reads2_files}\n'
            # )  # for debugging                
            CRTFsample_ = CRTFsample(name=name, genomeVersion=genomeVersion, reads1_files=reads1_files, reads2_files=reads2_files)
            if not CRTFsample_.source_files_check():
                continue
            CRTFsample_.get_input_reads_info()
            CRTFsample_.get_mapped_info()
            CRTFsample_.get_fragments_num()
            CRTFsample_.cutsites_files_check()
            CRTFsample_.get_peak_info()
            CRTFsample_.get_OCR_peak_info()
            mapping_csv.writerow(CRTFsample_.output_line())

            fragments_length_file = f'4_basic_QC/{name}_fragments_length.txt'
            if CRTFsample_.file_check(fragments_length_file, 'fragments length'):
                with open(fragments_length_file, 'r') as fhd:
                    next(fhd)
                    fL = [0 for i in range(301)]
                    for line in fhd:
                        line = line.strip().split()
                        try:
                            fL[int(line[0])] = int(line[1])
                        except IndexError:
                            pass
                    fL.insert(0, name)
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
        sys.stdout.write("User interrupts me! ;-) See you ^.^!")
        sys.exit(0)
