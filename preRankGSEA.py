#! /usr/bin/env python

import os, sys
from optparse import OptionParser, OptionGroup
from math import ceil, log10
import numpy as np
import pandas as pd
import time
import logging

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'Helvetica'
import matplotlib.pyplot as plt
import seaborn as sns

# ------------------------------------
# constants
# ------------------------------------
def logging_set(name):
    logging.basicConfig(level=10,
                format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                datefmt='%a, %d %b %Y %H:%M:%S',
                filename=name,
                filemode='w'
                )


# ------------------------------------
# Misc functions
# ------------------------------------
error   = logging.error    # function alias
warn    = logging.warning
debug   = logging.debug
info    = logging.info


# ------------------------------------
# Sub Functions
# ------------------------------------
def prepare_optparser():
    
    '''
    Prepare optparser object. New options will be added in thisfunction first.
    '''
    
    script_name = os.path.basename(sys.argv[0])
    usage = f'USAGE: {script_name} <-g geneSet> <-r rankedList> [-n name] [-s enrichmentStatistic]'
    description = 'pre-Ranked GSEA.'
    
    # option processor
    optparser = OptionParser(version=f'{script_name} 1.0', description=description, usage=usage, add_help_option=False)
    
    # basic setting
    optparser.add_option('-h','--help', action='help', help='Show this help message and exit.')
    optparser.add_option('-n', '--name',dest='name', type='string',
                         help='Name of this run. If not given, the time will be used. The output files will in the directory of this name under the running directory.')
    
    # input options
    group_input = OptionGroup(optparser, 'Input files options')
    group_input.add_option('-g', '--geneSet', dest='geneSet', type='string',
                              help='Interest gene set.')
    group_input.add_option('-r', '--rank', dest='rank', type='string',
                           help='Ranked gene list with given value.')
    optparser.add_option_group(group_input)
    
    # statistic options
    group_stats = OptionGroup(optparser, 'Statistic options')
    group_stats.add_option('-s', '--statistic', dest='statistic', type='float',
                           help='Weight for the genes with their value in rank. 0 means standard Kolmogovov-Smirnov statistic. 1 means weighting the genes by their value in rank.', default=1)
    # group_stats.added_option('-p', '--permutations', dest='permutations', type='int',
    #                          help='Number of permutations.', default=1000)
    optparser.add_option_group(group_stats)
    
    # plot options
    group_plot = OptionGroup(optparser, 'Plot options')
    group_plot.add_option('-p', '--plot', dest='plot', action='store_true',
                          help='Plot option. If ture, the program will generate ES enrichment plot.')
    group_plot.add_option('--xlim', dest='xlim',
                          help='X axis limits. --ylim 0,1')   
    group_plot.add_option('--ylim', dest='ylim',
                          help='Y axis limits for enrichment score. --ylim 0,1')    
    optparser.add_option_group(group_plot)
    
    return optparser


def opt_validate(optparser):
    
    '''Validate options from a OptParser object.

    Ret: Validated options object.
    '''
    (options,args) = optparser.parse_args()
    
    # input files must be given
    if not (options.geneSet and options.rank):
        error('Input geneSet and rank must be given!\n')
        optparser.print_help()
        sys.exit(1)
    
    if not os.path.isfile(options.geneSet):
        error('Can not find geneSet file: %s!\n' %options.geneSet)
        optparser.print_help()
        sys.exit(1)
    
    if not os.path.isfile(options.rank):
        error('Can not find rank file: %s!\n' %options.rank)
        optparser.print_help()
        sys.exit(1)
    
    if options.xlim:
        options.xlim = [float(value) for value in options.xlim.split(',')]
    if options.ylim:
        options.ylim = [float(value) for value in options.ylim.split(',')]
    
    # get name
    if not options.name:
        options.name='perRankGSEA_{}'.format(time.strftime('%Y_%b_%d_%H_%M_%S', time.localtime()))
    
    return options


def load_geneSet(geneSetFile):
    
    geneSet = []
    geneSetFhd = open(geneSetFile)
    for line in geneSetFhd:
        if line == 'Gene\n':
            continue
        geneSet.append(line.strip())
    geneSetFhd.close()
    
    return geneSet


def load_rank(rankFile):
    
    rank = []
    rankFhd = open(rankFile)
    for line in rankFhd:
        line = line.strip().split()
        rank.append((line[0], float(line[1])))
    rankFhd.close()
    rank.sort(key=lambda x:x[1], reverse=True)
    
    return rank


def enrichment_score(rank, geneSet, p):
    
    '''
    Calculate enrichment score for given rank & geneSet
    '''
    
    non_hit_gene = len(rank) - len(geneSet)
    Nr = sum([abs(rank[i][1])**p for i in range(len(rank)) if rank[i][0] in geneSet])
    
    pHit, pMiss, ES, hit = 0, 0, [0], []
    for i in range(len(rank)):
        if rank[i][0] in geneSet:
            pHit += abs(rank[i][1]) ** p / Nr
            hit.append(i)
        else:
            pMiss += 1.0 / non_hit_gene
        ES.append(pHit-pMiss)
    
    return (ES, hit)


def GSEA_plot(ES, hits, name):
    from matplotlib.offsetbox import AnchoredText
    with sns.axes_style('whitegrid'), sns.plotting_context('paper'):
        fig, ax = plt.subplots()
        ax.plot(np.arange(len(ES)),ES)
        maximum, minimum = max(ES), min(ES)
        ES_score = maximum if maximum+minimum>0 else minimum
        ylim_upper = maximum + (maximum - minimum) * .05
        ylim_lower = minimum - (maximum - minimum) * .05
        eventplot_height = (maximum - minimum) * .1
        ax.eventplot(hits, lineoffsets=ylim_lower-eventplot_height*.5,linelengths=eventplot_height)
        ax.set_xlim(0,len(ES)-1)
        ax.set_ylim((ylim_lower-eventplot_height,ylim_upper))
        ax.set_ylabel('Enrichment score')
        for ytick, yticklabel in zip(ax.get_yticks(),ax.get_yticklabels()):
            if ytick < ylim_lower:
                yticklabel.set_visible(False)
        ax.axhline(y=ylim_lower,xmin=0,xmax=len(ES)-1,color='k')
        at = AnchoredText(f'Enrichment score = {ES_score:.3f}', loc='upper right')
        ax.add_artist(at)
        plt.savefig(f'{name}.pdf')


# ------------------------------------
# Main Functions
# ------------------------------------

def main():

    # read the options and validate them
    options = opt_validate(prepare_optparser())

    # log setting
    if not os.path.isdir(options.name):
        os.mkdir(options.name)
    logging_set(f'{options.name}/{options.name}.log')

    # load data
    geneSet = load_geneSet(options.geneSet)
    rank = load_rank(options.rank)

    # calculate enrichment score
    (ES, hit) = enrichment_score(rank, geneSet, options.statistic)

    # output
    rFhd = open(options.name+'/ES.txt','w')
    rFhd.write('\n'.join([str(x) for x in ES])+'\n')
    rFhd.close()

    # plot
    if options.plot:
        GSEA_plot(ES, hit, f'{options.name}/{options.name}')


# ------------------------------------
# Program running
# ------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        warn('User interrupts me! ;-) See you ^.^!')
        sys.exit(0)

