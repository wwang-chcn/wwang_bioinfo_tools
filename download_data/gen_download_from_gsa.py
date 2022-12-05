#! /usr/bin/env python3

# ------------------------------------
# Load Modules
# ------------------------------------
import os
import sys
import optparse
from optparse import OptionParser
from typing import Iterable

import requests
from bs4 import BeautifulSoup

# ------------------------------------
# Constant
# ------------------------------------
root_url = 'https://download.cncb.ac.cn/'
data_root = 'gsa/'
ascp_server_base = 'aspera01@download.cncb.ac.cn:'


# ------------------------------------
# Sub Functions
# ------------------------------------
def accession_number_format_checker(accession_number: str) -> bool:
    if accession_number.startswith('CRA'):
        return accession_number[3:].isdigit()
    else:
        return False


def prepare_optparser() -> OptionParser:
    """
    Prepare optparser object. New options will be added in thisfunction first.
    """
    usage = "USAGE: %prog -a <accession number> -e <ascp path> -k <key file path> -o <output file name>"
    description = "gen_download_from_gsa -- Generate download scripts from GSA database."

    # option processor
    optparser = OptionParser(version="%prog 1.0",
                             description=description,
                             usage=usage,
                             add_help_option=True)
    # input options
    optparser.add_option(
        "-a",
        "--accession",
        dest="accession",
        type="string",
        action='append',
        help=
        "Required. Accession number. Multiple accession numbers can be input as -a accession_number_1 -aaccession_number_2."
    )
    optparser.add_option("-e",
                         "--execute",
                         dest="execute",
                         type="string",
                         help="Executeable ascp file path.")
    optparser.add_option("-k",
                         "--key",
                         dest="key",
                         type="string",
                         help="Ascp Key file path.")
    optparser.add_option("-o",
                         "--output",
                         dest="output",
                         type="string",
                         help="Output file name.")
    return optparser


def opt_validate(optparser: OptionParser) -> optparse.Values:
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """

    (options, _) = optparser.parse_args()

    # input accession numbers
    if options.accession is None:
        sys.stdout.write('Input accession numbers must be given!\n')
        optparser.print_help()
        sys.exit(1)
    for accession_number in options.accession:
        if not accession_number_format_checker(accession_number):
            sys.stdout.write(
                f'Input accession number ({accession_number}) is not valid.\n')
            optparser.print_help()
            sys.exit(1)

    # input executeable ascp
    if options.execute is None:
        sys.stdout.write('Input executeable ascp file must be given!\n')
        optparser.print_help()
        sys.exit(1)
    if not os.path.isfile(os.path.expanduser(options.execute)):
        sys.stdout.write(
            f'Can not find executeable ascp file: {options.execute}.\n')
        optparser.print_help()
        sys.exit(1)
    if not os.access(os.path.expanduser(options.execute), os.X_OK):
        sys.stdout.write(
            f'Input executeable ascp file ({options.execute}) is not executeable.\n'
        )
        optparser.print_help()
        sys.exit(1)

    # ascp key file
    if options.key is None:
        sys.stdout.write('Input ascp key file must be given!\n')
        optparser.print_help()
        sys.exit(1)
    if not os.path.isfile(os.path.expanduser(options.key)):
        sys.stdout.write(f'Can not find ascp key file: {options.key}.\n')
        optparser.print_help()
        sys.exit(1)

    # output file
    if options.output is None:
        sys.stdout.write('Output file must be given!\n')
        optparser.print_help()
        sys.exit(1)
    if os.path.isfile(os.path.expanduser(options.output)):
        sys.stdout.write(
            f'Warning! Output file ({options.output}) has already exist.\n')

    return options


def gen_run_accession(cra_accessions: Iterable[str]) -> Iterable[str]:
    """Generate run accession based on given cra accession.
Parameters
----------
cra_accessions : Iterable[str]
                 CRA accessions.

Returns
----------
run_accessions : Iterable[str]
                 Iterable Run accessions.
"""
    for cra_accession in cra_accessions:
        cra_url = root_url + data_root + f'{cra_accession}/'
        for node in BeautifulSoup(requests.get(cra_url).text,
                                  'html.parser').find_all('a'):
            path = node.get('href')
            if path.startswith('CRR'):
                yield f'{cra_accession}/' + path


def gen_file_url(run_accessions: Iterable[str]) -> Iterable[str]:
    """Generate sequencing file url based on given Run accession.
Parameters
----------
run_accessions : Iterable[str]
                 Iterable Run accessions.

Returns
----------
file_urls : Iterable[str]
            Iterable sequencing file urls.
"""
    for run_accession in run_accessions:
        run_url = root_url + data_root + run_accession
        for node in BeautifulSoup(requests.get(run_url).text,
                                  'html.parser').find_all('a'):
            path = node.get('href')
            if path.endswith('.gz'):
                yield ascp_server_base + data_root + run_url + path


def gen_download_scrip(ascp_file: str, key_file: str,
                       file_urls: Iterable[str]) -> Iterable[str]:
    """Generate download script based on given file url.
Parameters
----------
file_urls : Iterable[str]
            Iterable sequencing file urls.

Returns
----------
download_scirpts : Iterable[str]
                   Iterable sequencing file download script.
"""
    for file_url in file_urls:
        yield f'{ascp_file} -QT -l 300m -P33001 -i {key_file} {file_url} {file_url.rsplit("/", 1)[1]}'


# ------------------------------
# Main Funcitons
# ------------------------------
def main():
    """Main function"""

    # read the options and validate them
    options = opt_validate(prepare_optparser())

    run_accessions = gen_run_accession(options.accession)
    file_urls = gen_file_url(run_accessions)
    download_scirpts = gen_download_scrip(options.execute, options.key,
                                          file_urls)
    with open(options.output, 'w') as fhd:
        for download_scirpt in download_scirpts:
            fhd.write(download_scirpt + '\n')


# ------------------------------
# Program Running
# ------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stdout.write('User interrupts me! ;-) See you!\n')
        sys.exit(0)
