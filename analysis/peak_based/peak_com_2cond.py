import os
import subprocess
import sys
from copy import deepcopy
from optparse import OptionParser, Values
from typing import Dict, List

import matplotlib as mpl
import numpy as np
import pandas as pd

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'Arial'

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2

# ------------------------------------
# Constant
# ------------------------------------
FILE_EXT: Dict[str, List[str]] = {
    'narrowPeak': ['_peaks.narrowPeak'],
    'bigWig': ['.bw', '.bigWig'],
}

ENV_PREFIX = '''source ~/.bash_profile
conda activate hs
'''

BASE_PATH = os.getcwd()

SELF_DIR = os.path.abspath(os.path.dirname(__file__))


# ------------------------------------
# functions
# ------------------------------------
def prepare_optparser():
    '''\
    Prepare optparser object. New options will be added in thisfunction first.
    '''
    script_name = os.path.basename(sys.argv[0])
    usage = f'USAGE: {script_name} <-n NAME> <-p PATH> <--c1 CONDITION1> <--c2 CONDITION2> <--c1-name CONDITION1-NAME> <--c2-name CONDITION2-NAME> <-g GENOME> [-o OUT=.] [-q Q_CUTOFF=6] [-f F_CUTOFF=6] [-high HIGH=10] [-low LOW=2]'
    description = 'peak_com_2cond -- merge peaks from replicates'

    # option processor
    optparser = OptionParser(version=f'{script_name} 0.1', description=description, usage=usage, add_help_option=False)

    # basic setting
    optparser.add_option('-h', '--help', action='help', help='Show this help message and exit.')
    optparser.add_option('-n', '--name', dest='name', type='string', help='Name for this compare.')
    optparser.add_option('-p',
                         '--path',
                         dest='path',
                         type='string',
                         default='.',
                         help='Path where store the peak and bigWigfile.')
    optparser.add_option('--c1',
                         '--condition1',
                         dest='condition1',
                         type='string',
                         help='Name base of condition 1 peak and bigWig file.')
    optparser.add_option('--c2',
                         '--condition2',
                         dest='condition2',
                         type='string',
                         help='Name base of condition 2 peak and bigWig file.')
    optparser.add_option('--c1-name',
                         '--condition1-name',
                         dest='condition1_name',
                         type='string',
                         help='Name for condition 1.')
    optparser.add_option('--c2-name',
                         '--condition2-name',
                         dest='condition2_name',
                         type='string',
                         help='Name for condition 2.')
    optparser.add_option('-g', '--genome', dest='genome', type='string', help='Genome name.')
    optparser.add_option('-o',
                         '--out',
                         dest='out',
                         type='string',
                         default='.',
                         help='Path for output files. A subfolder with the name will be created in the path.')
    optparser.add_option('-q', '--q-cutoff', dest='q_cutoff', type='float', help='Qvalue cutoff.', default='6')
    optparser.add_option('-f', '--f-cutoff', dest='f_cutoff', type='float', help='Fold cutoff.', default='6')
    optparser.add_option('--high', dest='high', type='float', help='High threshold for signal.', default='10')
    optparser.add_option('--low', dest='low', type='float', help='Low threshold for signal.', default='2')

    return optparser


def opt_validate(optparser) -> Values:
    '''Validate options from a OptParser object.
    Ret: Validated options object.
    '''
    (options, args) = optparser.parse_args()

    # input name must be given
    if not options.name:
        sys.stdout.write('Input name must be given!\n')
        optparser.print_help()
        sys.exit(1)

    # input path must be given
    if not options.path:
        sys.stdout.write('Input path must be given!\n')
        optparser.print_help()
        sys.exit(1)
    # check the path
    if not os.path.isdir(options.path):
        sys.stdout.write(f'Input path ({options.path}) must be a directory!\n')
        optparser.print_help()
        sys.exit(1)

    # input condition1 and condition2 must be given
    if not (options.condition1 and options.condition2):
        sys.stdout.write('Condition 1 and 2 must be given!\n')
        optparser.print_help()
        sys.exit(1)
    # check the existence of peak and bigWig files
    file_type = ['narrowPeak', 'bigWig']
    for cond in [options.condition1, options.condition2]:
        for ft in file_type:
            exts = FILE_EXT[ft]
            for ext in exts:
                if os.path.isfile(f'{options.path}/{cond}{ext}'):
                    break
            else:
                sys.stdout.write(f'Can not find the {ft} file: {options.path}/{cond}{exts[0]}.\n')
                optparser.print_help()
                sys.exit(1)

    # input condition1_name and condition2_name must be given
    if not (options.condition1_name and options.condition2_name):
        sys.stdout.write('Condition 1 and 2 name must be given!\n')
        optparser.print_help()
        sys.exit(1)

    return options


def load_input_files(path: str, condition1: str, condition2: str, condition1_name: str, condition2_name: str) -> None:
    """
    Create symlink files for the peak and bigWig files for the two conditions.
    The function supporsed to be called under the output folder.
    Peak and bigWigs files should be checked before calling this function.
    """
    file_type = ['narrowPeak', 'bigWig']
    for cond, cond_name in zip([condition1, condition2], [condition1_name, condition2_name]):
        for ft in file_type:
            exts = FILE_EXT[ft]
            for ext in exts:
                src_file = os.path.join(BASE_PATH, path, f'{cond}{ext}')
                print(src_file)
                if os.path.isfile(src_file):
                    print('src_file exists')
                    dst_file = f'{cond_name}{ext}'
                    if not os.path.islink(dst_file):
                        print(f'create symlink: {dst_file}')
                        os.symlink(src=os.path.relpath(src_file), dst=dst_file)
                    break


def merge_peaks(name: str, condition1_name: str, condition2_name: str, q_cutoff: float, f_cutoff: float,
                genome: str) -> None:
    """
    Merge peaks from replicates.
    """
    cmd = ENV_PREFIX
    peak_file_option = ' '.join([f'-p {condition1_name}_peaks.narrowPeak', f'-p {condition2_name}_peaks.narrowPeak'])
    bigwig_file_option = ' '.join([f'-b {condition1_name}.bw', f'-b {condition2_name}.bw'])
    genome_size_file_option = os.path.expanduser(f'-g ~/source/bySpecies/{genome}/{genome}.chrom.sizes')
    cmd += f'python {SELF_DIR}/merge_peak.py -n {name} {peak_file_option} -q {q_cutoff} -f {f_cutoff} {bigwig_file_option} {genome_size_file_option}\n'
    cmd += f'intersectBed -c -a {name}_summits.bed -b {condition1_name}_peaks.narrowPeak | intersectBed -c -a - -b {condition2_name}_peaks.narrowPeak > {name}_merged_peak_overlap.bed\n'
    subprocess.call(cmd, shell=True)


def generate_summit_100bp(name: str, condition_1_name: str, condition_2_name: str) -> pd.DataFrame:
    """
    Generate 100bp summit file for signal capture.
    """

    peak_overlap = pd.read_csv(f"{name}_merged_peak_overlap.bed", sep='\t', header=None)
    peak_overlap.columns = ['chr', 'summit', 'summit+1', 'name', 'score', 'strand', condition_1_name, condition_2_name]

    peak_overlap_100bp = peak_overlap[['chr', 'summit', 'name']]
    peak_overlap_100bp['start'] = peak_overlap_100bp['summit'] - 50
    peak_overlap_100bp['end'] = peak_overlap_100bp['summit'] + 50
    peak_overlap_100bp = peak_overlap_100bp[['chr', 'start', 'end', 'name']]
    peak_overlap_100bp.to_csv(f"{name}_summits_100bp.bed", sep='\t', index=None, header=None)

    return peak_overlap


def capture_signal_over_merged_peaks(name: str, condition_1_name: str, condition_2_name: str) -> None:
    """
    Capture signal over merged peaks."""
    cmd = ENV_PREFIX
    cmd += f'bigWigAverageOverBed {condition_1_name}.bw {name}_summits_100bp.bed {name}_summits_100bp_{condition_1_name}.tsv &\n'
    cmd += f'bigWigAverageOverBed {condition_2_name}.bw {name}_summits_100bp.bed {name}_summits_100bp_{condition_2_name}.tsv &\n'
    cmd += 'wait\n'
    subprocess.call(cmd, shell=True)


def diff_peak_status(name: str, condition_1_name: str, condition_2_name: str, peak_over_df: pd.DataFrame,
                     high_threshold: float, low_threshold: float) -> pd.DataFrame:
    """
    Diff peak between two conditions.
    """

    def status_gen():
        for _, row in diff_signal_df.iterrows():
            if row['ov_status'] == 'both':
                if row[condition_1_name] >= high_threshold and row[condition_2_name] >= high_threshold:
                    yield 'both'
                else:
                    yield 'other'
            elif row['ov_status'] == f'{condition_1_name}-specific':
                if row[condition_1_name] >= high_threshold and row[condition_2_name] < low_threshold:
                    yield f'{condition_1_name}-specific'
                else:
                    yield 'other'
            elif row['ov_status'] == f'{condition_2_name}-specific':
                if row[condition_2_name] >= high_threshold and row[condition_1_name] < low_threshold:
                    yield f'{condition_2_name}-specific'
                else:
                    yield 'other'

    signal_df_dict: Dict[str, pd.DataFrame] = {
        condition_1_name:
        pd.read_csv(f"{name}_summits_100bp_{condition_1_name}.tsv",
                    sep='\t',
                    names=['name', 'size', 'covered', 'sum', 'mean0', 'mean']),
        condition_2_name:
        pd.read_csv(f"{name}_summits_100bp_{condition_2_name}.tsv",
                    sep='\t',
                    names=['name', 'size', 'covered', 'sum', 'mean0', 'mean'])
    }
    diff_signal_df = signal_df_dict[condition_1_name][['name', 'mean0']]
    diff_signal_df.columns = ['name', condition_1_name]
    diff_signal_df[condition_2_name] = signal_df_dict[condition_2_name]['mean0']
    diff_signal_df['ov_status'] = [
        'both' if x and y else (f'{condition_1_name}-specific' if x and not y else f'{condition_2_name}-specific')
        for x, y in zip(peak_over_df[condition_1_name], peak_over_df[condition_2_name])
    ]
    diff_signal_df['peak_status'] = list(status_gen())

    return diff_signal_df


def visualize_diff_peak_status(diff_signal_df: pd.DataFrame, name: str, condition_1_name: str,
                               condition_2_name: str) -> None:
    with sns.axes_style('white', rc={
            'xtick.bottom': True,
            'ytick.left': True
    }), sns.plotting_context('paper',
                             rc={
                                 'axes.titlesize': 8,
                                 'axes.labelsize': 8,
                                 'xtick.labelsize': 6,
                                 'ytick.labelsize': 6,
                                 'legend.fontsize': 6
                             }):

        # Scatter plot

        fig, ax = plt.subplots(figsize=(4, 4))
        sns.scatterplot(data=diff_signal_df,
                        x=condition_1_name,
                        y=condition_2_name,
                        hue='peak_status',
                        s=10,
                        edgecolor=None,
                        ax=ax)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(f'{condition_1_name} signal')
        ax.set_ylabel(f'{condition_2_name} signal')
        fig.tight_layout()
        fig.savefig(f'Scatter_{name}_diff_peak_status.pdf')

        # Ven diagram
        condition_1_specific_num = sum(diff_signal_df['peak_status'] == f'{condition_1_name}-specific')
        condition_2_specific_num = sum(diff_signal_df['peak_status'] == f'{condition_2_name}-specific')
        both_num = sum(diff_signal_df['peak_status'] == f'both')

        fig, ax = plt.subplots(figsize=(4, 4))
        venn2(subsets=(condition_1_specific_num, condition_2_specific_num, both_num),
              set_labels=(condition_1_name, condition_2_name),
              ax=ax)
        fig.tight_layout()
        fig.savefig(f'Venn_{name}_peak_num.pdf', transparent=True)


def avg_signal_around_sites(bed_file: str,
                            bigwig_files: List[str],
                            name: str,
                            labels: List[str],
                            resolution: int,
                            span: int,
                            normalization: bool = True,
                            strand: bool = True,
                            bw_scan: bool = True) -> Dict[str, pd.DataFrame]:
    import os
    import subprocess

    # check strand parameters
    bed_fields = int(subprocess.check_output(f'''head -1 {bed_file} | awk '{{print NF}}' ''', shell=True).decode())
    strand = strand & (bed_fields >= 6)

    # useful funcs
    def get_bigwig_mean(bw_file):
        import subprocess
        cmd = ENV_PREFIX
        cmd += f'''bigWigInfo {bw_file} | grep mean'''
        fold = subprocess.check_output(args=cmd, shell=True).decode()
        fold = float(fold.split()[1])
        return fold

    xticks = np.arange(-span, span + resolution, resolution)
    capture_points = int(2 * span / resolution + 1)

    # capture regions file prepare
    if strand:
        cmd = f'''awk '{{mid=($2+$3)/2; printf $1"\\t%d\\t%d\\t"$4"\\t"$5"\\t"$6"\\n", mid-{span}-{resolution >> 1}, mid+{span}+{resolution >> 1}}}' {bed_file} > capture_regions.bed'''
    else:
        cmd = f'''awk '{{mid=($2+$3)/2; printf $1"\\t%d\\t%d\\t"$4"\\n", mid-{span}-{resolution >> 1}, mid+{span}+{resolution >> 1}}}' {bed_file} > capture_regions.bed'''
    print(f'cmd: {cmd}')
    print(subprocess.check_output(cmd, shell=True).decode(), end='')

    # capture signal
    capture_regions_file = 'capture_regions.bed'
    bw_scan_cmd = ENV_PREFIX
    bw_scan_cmd += f'''python {SELF_DIR}/../../getBigWigValue.py -n {name} -p 4 -s {capture_points} -m speed -b {capture_regions_file}'''
    for bigwig_file in bigwig_files:
        bw_scan_cmd += f' -w {bigwig_file}'
    if strand:
        bw_scan_cmd += ' --strand'
    if bw_scan:
        print(f'bw_scan_cmd: {bw_scan_cmd}')
        print(subprocess.check_output(bw_scan_cmd, shell=True).decode(), end='')
    os.remove('capture_regions.bed')

    # output
    signal = {}
    for index, (bigwig_file, label) in enumerate(zip(bigwig_files, labels)):
        mean = get_bigwig_mean(bigwig_file) if normalization else 1
        signal[label] = pd.read_csv(
            f'signal_{name}_siteprof{index+1}.gz', sep='\t', header=None, index_col=None, names=xticks) / mean

    return signal


def capture_signal_around_sites(name: str,
                                condition_1_name: str,
                                condition_2_name: str,
                                resolution: int = 10,
                                span: int = 1000) -> Dict[str, pd.DataFrame]:
    """
    Capture signal around sites.
    """
    signal_df_dict = avg_signal_around_sites(
        bed_file=f"{name}_summits.bed",
        bigwig_files=[f"{condition_1_name}.bw", f"{condition_2_name}.bw"],
        name=name,
        labels=[condition_1_name, condition_2_name],
        resolution=resolution,
        span=span,
        strand=False,
    )
    return signal_df_dict


def avg_signal_around_sites_plot(signal,
                                 group=None,
                                 group_order=None,
                                 xlabel='Fold',
                                 ylabel='Distance from center (bp)',
                                 common_ylim=True,
                                 output_file='lineplot_group.pdf'):
    signal_ = deepcopy(signal)
    if group is None:
        avg_signal = pd.DataFrame([])
        for label, signal_matrix in signal_.items():
            avg_signal_ = signal_matrix.mean()
            for pos in avg_signal_.index:
                avg_signal = pd.concat(
                    [avg_signal,
                     pd.DataFrame({
                         'pos': [pos],
                         'value': [avg_signal_[pos]],
                         'sample': [label]
                     })])
        with sns.axes_style('white', rc={
                'xtick.bottom': True,
                'ytick.left': True
        }), sns.plotting_context('paper',
                                 rc={
                                     'axes.titlesize': 14,
                                     'axes.labelsize': 12,
                                     'xtick.labelsize': 10,
                                     'ytick.labelsize': 10,
                                     'legend.fontsize': 10
                                 }):
            fig, ax = plt.subplots()
            sns.lineplot(x='pos', y='value', hue='sample', data=avg_signal, ax=ax)
            ax.legend()
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(label)
            fig.tight_layout()
            fig.savefig(output_file, transparent=True)

    else:
        group_order = list(set(group)) if group_order is None else group_order
        avg_signal = pd.DataFrame([])
        for label, signal_matrix in signal_.items():
            signal_matrix['group'] = group
            for index, row in signal_matrix.groupby('group').mean().iterrows():
                for pos in row.index:
                    avg_signal = pd.concat([
                        avg_signal,
                        pd.DataFrame({
                            'group': [index],
                            'pos': [pos],
                            'value': [row[pos]],
                            'sample': [label]
                        })
                    ])
        with sns.axes_style('white', rc={
                'xtick.bottom': True,
                'ytick.left': True
        }), sns.plotting_context('paper',
                                 rc={
                                     'axes.titlesize': 14,
                                     'axes.labelsize': 12,
                                     'xtick.labelsize': 10,
                                     'ytick.labelsize': 10,
                                     'legend.fontsize': 10
                                 }):
            fig, ax = plt.subplots(len(signal_), 1, figsize=(6.4, 4.8 * 0.6 * len(signal_)))
            for index, label in enumerate(signal_.keys()):
                sns.lineplot(x='pos',
                             y='value',
                             hue='group',
                             hue_order=group_order,
                             data=avg_signal.loc[avg_signal['sample'] == label, :],
                             ax=ax[index])
                ax[index].legend()
                ax[index].set_xlabel(xlabel)
                ax[index].set_ylabel(ylabel)
                ax[index].set_title(label)
            if common_ylim:
                ymin, ymax = [], []
                for index in range(len(signal_)):
                    ylim = ax[index].get_ylim()
                    ymin.append(ylim[0])
                    ymax.append(ylim[1])
                ymin, ymax = min(ymin), max(ymax)
                for index in range(len(signal_)):
                    ax[index].set_ylim((ymin, ymax))
            fig.tight_layout()
            fig.savefig(output_file, transparent=True)


def avg_signal_around_sites_plot_group(signal_df_dict: Dict[str, pd.DataFrame], diff_signal_df: pd.DataFrame, name: str,
                                       condition_1_name: str, condition_2_name: str) -> None:
    avg_signal_around_sites_plot(signal_df_dict,
                                 group=diff_signal_df['peak_status'],
                                 group_order=['both', f'{condition_1_name}-specific', f'{condition_2_name}-specific'],
                                 xlabel='Distance from motif center (bp)',
                                 ylabel='Fold of ATAC-seq signal',
                                 output_file=f'Lineplot_{name}_signal_around_peak_summits.pdf')


def signal_heatmap(signal,
                   span,
                   output_file='heatmap.pdf',
                   ranking_signal=None,
                   group=None,
                   group_order=None,
                   group_separate=False,
                   group_splitter=False,
                   vmin=0,
                   vmax=30,
                   cmap=None,
                   figsize=None,
                   height=None,
                   xlabel='Distance from peak summits (bp)'):
    # load modules
    from itertools import cycle, islice

    signal_num = len(signal)
    # handle parameters
    # #vmin&vmax
    if isinstance(vmin, (float, int)):
        vmin = {k: vmin for k in signal}
    if isinstance(vmax, (float, int)):
        vmax = {k: vmax for k in signal}
    # #cmap
    if cmap is None:
        cmap = ['RdBu_r'] * signal_num
    elif isinstance(cmap, list):
        cmap = list(islice(cycle(cmap), signal_num))
    else:
        cmap = [cmap] * signal_num
    # #figsize
    if figsize is None:
        figsize = (3 * signal_num, 8)
    # #order: (ranking_signal & group)
    if group is None and ranking_signal is None:  # no order
        order = None
    else:
        order_df = pd.DataFrame()
        if group is not None:
            if group_order is None:
                group_order = list(set(group))
            group_index = {group_name: i for i, group_name in enumerate(group_order)}
            order_df['group'] = [group_index[g] if g in group_index else len(group_index) for g in group]
        if ranking_signal is not None:
            order_df['ranking_signal'] = ranking_signal
            order_df.sort_values(by='ranking_signal', ascending=False, inplace=True)
        if group is not None:
            group_nums = [0] + list(np.cumsum([group.count(group_name) for group_name in group_order]))
            order_df.sort_values(by='group', inplace=True)
        order = order_df.index.to_numpy()

    # xticks
    if span >= 1000:
        xticks = np.arange(-span, span + 1000, 1000)
        xticklabels = [f'{x/1000: .0f}' for x in xticks]
        xticklabels[-1] = f'{xticklabels[-1]} kb'
    elif span >= 100:
        xticks = np.arange(-span, span + 100, 100)
        xticklabels = [f'{x: d}' for x in xticks]
    else:
        xticks = np.arange(-span, span + 10, 10)
        xticklabels = [f'{x: d}' for x in xticks]
    if height is None:
        height = 3000

    with sns.axes_style('white'), sns.plotting_context('notebook',
                                                       rc={
                                                           'axes.titlesize': 8,
                                                           'axes.labelsize': 8,
                                                           'xtick.labelsize': 6,
                                                           'ytick.labelsize': 6,
                                                       }):
        if group is not None and group_separate:
            fig, axs = plt.subplots(len(group_order), signal_num, figsize=figsize)
            for index, (label, signal_matrix) in enumerate(signal.items()):
                for group_index, group_name in enumerate(group_order):
                    signal_matrix_ = signal_matrix.values[order[group_nums[group_index]:group_nums[group_index + 1]], :]
                    c = axs[group_index, index].imshow(signal_matrix_,
                                                       cmap=cmap[index],
                                                       extent=[-span, span, height, 0],
                                                       interpolation='nearest',
                                                       origin='upper',
                                                       vmin=vmin[label],
                                                       vmax=vmax[label])
                    if group_index == len(group_order) - 1:
                        axs[group_index, index].set_xticks(xticks)
                        axs[group_index, index].set_xticklabels(xticklabels)
                        axs[group_index, index].set_xlabel(xlabel)
                    if index == 0:
                        axs[group_index, index].set_ylabel(group_order[group_index])
                        axs[group_index, index].set_yticks([])
                        axs[group_index, index].set_yticklabels([])
                    if group_index == 0:
                        axs[group_index, index].set_title(label)
                    # fig.colorbar(c, ax=axs[group_index,index], orientation='horizontal')
            fig.savefig(output_file, transparent=True)
        else:
            fig, axs = plt.subplots(1, signal_num, figsize=figsize)
            for index, (label, signal_matrix) in enumerate(signal.items()):
                signal_matrix = signal_matrix.values[order, :]
                c = axs[index].imshow(signal_matrix,
                                      cmap=cmap[index],
                                      extent=[-span, span, height, 0],
                                      interpolation='nearest',
                                      origin='upper',
                                      vmin=vmin[label],
                                      vmax=vmax[label])
                axs[index].set_xticks(xticks)
                axs[index].set_xticklabels(xticklabels)
                axs[index].set_yticks([])
                axs[index].set_yticklabels([])
                axs[index].set_xlabel(xlabel)
                axs[index].set_title(label)
                if group is not None and group_splitter:
                    axs[index].hlines(y=[x + 0.5 for x in group_nums[1:-1]], xmin=-span, xmax=span, ls='--', lw=2)
                # fig.colorbar(c, ax=axs[index], orientation='horizontal')
            fig.tight_layout()
            fig.savefig(output_file, transparent=True)


def avg_signal_around_sites_heatmap(signal_df_dict: Dict[str, pd.DataFrame], diff_signal_df: pd.DataFrame, name: str,
                                    condition_1_name: str, condition_2_name: str) -> None:

    output_file = f'Heatmap_{name}_signal_around_peak_summits.pdf'
    group = diff_signal_df['peak_status'].tolist()
    group_order = ['both', f'{condition_1_name}-specific', f'{condition_2_name}-specific']
    colors = ['#FFFFFF', '#CD2626']
    cmap1 = mpl.colors.LinearSegmentedColormap.from_list("mycmap", colors)
    #
    signal_heatmap(
        signal=signal_df_dict,
        span=list(signal_df_dict.values())[0].columns[-1],
        output_file=output_file,
        ranking_signal=np.arange(diff_signal_df.shape[0]),
        group=group,
        group_order=group_order,
        group_separate=True,
        group_splitter=False,
        vmin=0,
        vmax=10,
        cmap=cmap1,
        figsize=None,
        height=None,
        xlabel='Distance from peak summits (bp)',
    )


# ------------------------------------
# Main function
# ------------------------------------
def main() -> None:

    # read the options and validate them
    options: Values = opt_validate(prepare_optparser())

    # create the output folder
    out_dir = os.path.join(options.out, options.name)
    os.makedirs(out_dir, exist_ok=True)
    os.chdir(out_dir)

    # create symlink for the peak and bigWig files
    load_input_files(path=options.path,
                     condition1=options.condition1,
                     condition2=options.condition2,
                     condition1_name=options.condition1_name,
                     condition2_name=options.condition2_name)

    # merge peaks
    merge_peaks(name=options.name,
                condition1_name=options.condition1_name,
                condition2_name=options.condition2_name,
                q_cutoff=options.q_cutoff,
                f_cutoff=options.f_cutoff,
                genome=options.genome)

    # generate 100bp summit file
    peak_overlap_df = generate_summit_100bp(name=options.name,
                                            condition_1_name=options.condition1_name,
                                            condition_2_name=options.condition2_name)
    peak_overlap_df.to_csv(f"{options.name}_merged_peak_overlap.bed", sep='\t', index=None, header=None)

    # capture signal over merged peaks
    capture_signal_over_merged_peaks(name=options.name,
                                     condition_1_name=options.condition1_name,
                                     condition_2_name=options.condition2_name)

    # diff peak status
    diff_signal_df = diff_peak_status(name=options.name,
                                      condition_1_name=options.condition1_name,
                                      condition_2_name=options.condition2_name,
                                      peak_over_df=peak_overlap_df,
                                      high_threshold=options.high,
                                      low_threshold=options.low)
    diff_signal_df.to_csv(f"{options.name}_diff_peak_status.tsv", sep='\t', index=None)

    # Visualize the diff peak status
    visualize_diff_peak_status(diff_signal_df=diff_signal_df,
                               name=options.name,
                               condition_1_name=options.condition1_name,
                               condition_2_name=options.condition2_name)

    # Capture signal around sites
    signal_df_dict = capture_signal_around_sites(name=options.name,
                                                 condition_1_name=options.condition1_name,
                                                 condition_2_name=options.condition2_name)

    # Plot the signal around sites
    avg_signal_around_sites_plot_group(signal_df_dict=signal_df_dict,
                                       diff_signal_df=diff_signal_df,
                                       name=options.name,
                                       condition_1_name=options.condition1_name,
                                       condition_2_name=options.condition2_name)
    avg_signal_around_sites_heatmap(signal_df_dict=signal_df_dict,
                                    diff_signal_df=diff_signal_df,
                                    name=options.name,
                                    condition_1_name=options.condition1_name,
                                    condition_2_name=options.condition2_name)


# ------------------------------------
# Program Running
# ------------------------------------

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stdout.write('User interrupts me! ;-) See you!\n')
        sys.exit(0)
