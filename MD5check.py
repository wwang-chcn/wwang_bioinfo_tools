#! /usr/bin/env python3

# Nov-1-2018

# MD5 auto check on server

import subprocess
import sys
import threading
from distutils.spawn import find_executable
from optparse import OptionParser


def prepare_optparser():
    '''\
    Prepare optparser object. New options will be added in thisfunction first.
    '''
    usage = 'USAGE: %prog -p [processor]'
    description = 'bigWigSignalCapture -- Capture signal form bigWigFiles.'

    # option processor
    optparser = OptionParser(version='% prog 0.1',
                             description=description,
                             usage=usage,
                             add_help_option=True)

    # basic setting
    optparser.add_option('-p', '--process',dest='process', type='int',\
                         help='''Number of subprocess. Deafault number is 4.
Please set to 1 when checking data on low speed disk like external hard drive.''',default=4)
    return optparser


def plotform_check():
    operation_system = sys.platform
    if operation_system == "darwin":
        checksumCMD = find_executable("md5")
    elif operation_system.startswith("linux"):
        checksumCMD = find_executable("md5sum")
    else:
        sys.stdout.write(
            "Do not support operation system except LINUX or macOS, exit!\n")
        sys.exit(2)
    if not checksumCMD:
        sys.stdout.write("Can not find executable md5/md5sum, exit!\n")
        sys.exit(3)
    return checksumCMD


def parse_md5_result(result):
    global lock
    line = result.strip().split()
    name, md5 = '', ''
    if len(line) == 2:
        md5, name = line
    elif len(line) == 4:
        _, name, _, md5 = line
        name = name[1:-1]
    else:
        lock.acquire()
        sys.stdout.write('Undetermined md5checkoutput format: {}.\n'.format(
            '\t'.join(line)))
        lock.release()
    return (name.split('/')[-1], md5)


def md5_compare(name, md5, sampleMD5):
    lock.acquire()
    try:
        if md5 == sampleMD5[name]:
            sys.stdout.write('Sample: {} is valid.\n'.format(name))
        else:
            sys.stdout.write(
                'Sample: {} is not valid. Expect md5: {}. Observed md5: {}.\n'.
                format(name, sampleMD5[name], md5))
    except KeyError:
        sys.stdout.write(
            'Sample: {} is not valid. Can not find expect md5 information. The observed md5 is: {}.\n'
            .format(name, md5))
    lock.release()


def get_md5():
    sp = subprocess.Popen('{command} {args}'.format(
        command='cat', args='*md5*txt *md5 MD5* */*md5*txt */*md5 */MD5*'),
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          shell=True,
                          encoding='utf8')
    out, err = sp.communicate()
    sampleMD5 = dict(map(parse_md5_result, out.strip().split("\n")))
    return (sampleMD5)


def check(checksumCMD):
    global lock
    global failed_samples
    global index
    global fastq_files
    global sampleMD5
    while True:
        lock.acquire()
        if index < len(fastq_files):
            index += 1
            lock.release()
            sp = subprocess.Popen('{command} {args}'.format(
                command=checksumCMD, args=fastq_files[index - 1]),
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  shell=True,
                                  encoding='utf8')
            out, err = sp.communicate()
            name, md5 = parse_md5_result(out)
            md5_compare(name, md5, sampleMD5)
        else:
            lock.release()
            break


(options, args) = prepare_optparser().parse_args()
lock = threading.Lock()
checksumCMD = plotform_check()
sampleMD5 = get_md5()
sp = subprocess.Popen(" ".join(
    ["ls", "*/*fastq.gz", "*fastq.gz", "*/*fq.gz", "*fq.gz"]),
                      stdout=subprocess.PIPE,
                      stderr=subprocess.PIPE,
                      shell=True,
                      encoding='utf8')
out, err = sp.communicate()
fastq_files = out.strip().split()
index = 0
failed_samples = {}

for i in range(options.process):
    new_thread = threading.Thread(target=check, args=(checksumCMD, ))
    new_thread.start()

##def main():
##
##    checkMD5()
##
### ------------------------------------
### Program running
### ------------------------------------
##if __name__ == '__main__':
##    try:
##        main()
##    except KeyboardInterrupt:
##        warn("User interrupts me! ;-) See you ^.^!")
##        sys.exit(0)
