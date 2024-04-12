import pandas as pd
import pybedtools

import argparse
import os
import sys
import time
import datetime

import multiprocessing

from gtfparse import read_gtf
import warnings
import logging

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import seaborn as sns

from scipy.stats import mannwhitneyu
from statsmodels.stats import multitest


def parse_args():
    parser = argparse.ArgumentParser(description='''Characterize the depletion rank scores of genomic fragment class (genomic_fragment_class.txt) from PIN.py.''')
    
    parser.add_argument('-g', '--geno-frag-class', required=True, help='''Genomic fragment class file from PIN.py (genomic_fragment_class.txt).''')
    parser.add_argument('-dr', '--dr-ref', required=True, default=os.path.join(os.path.dirname(__file__), '../lib/reference_files/DR.gor'), help='''DR score reference file (DR.gor).''')
    parser.add_argument('-ga', '--geno_anno', required=True, default=os.path.join(os.path.dirname(__file__), '../lib/reference_files/gencode.v43.annotation.gtf.gz'), help='''GENCODE Genome annotation file (e.g. gencode.v43.annotation.gtf.gz).''')
    parser.add_argument('-o', '--output-dir', default=os.path.join(os.getcwd(), 'results'), help='''Directory to write results.''')
    parser.add_argument('-p', '--num-processes', type=int, default=int(multiprocessing.cpu_count()/2), help='Number of CPUs to use (default: half the number of CPUs).')

    return parser.parse_args()

class Logger:
    def __init__(self, output_dir):
        self.console = sys.stdout
        self.log = open(os.path.join(output_dir, 'depletion_rank.log'), 'w')

        time = datetime.datetime.now()
        self.write((f'{time.strftime("%A")},'
                    f' {time.strftime("%b")}' 
                    f' {time.strftime("%d")},'
                    f' {time.strftime("%Y")}' 
                    f' {time.strftime("%I")}:'
                    f'{time.strftime("%M")}'
                    f' {time.strftime("%p")}'))
        
    def write(self, message):
        self.console.write(message+'\n')
        self.log.write(message+'\n')
        self.log.flush()

    def close(self):
        self.log.close()

def process_input(dr_ref, geno_frag_class, geno_anno):
    dr = pd.read_csv(dr_ref, sep='\t')
    dr = dr.rename(columns={'Chr':'DR_chr', 'Fromx':'DR_start', 'To':'DR_end', 'rank':'DR_score'})
    dr['DR_start'] = dr['DR_start'].sub(1)

    geno_frag_class = pd.read_csv(geno_frag_class, sep='\t')

    gene_anno = process_genome_annotation(geno_anno)
    gene_anno_exon = gene_anno[gene_anno['feature']=='exon'][['seqname', 'start', 'end']].drop_duplicates()

    return dr, geno_frag_class, gene_anno_exon

def process_genome_annotation(geno_anno):
    # Suppress FutureWarnings
    warnings.simplefilter(action='ignore', category=FutureWarning)

    # Configure the root logger to only display warnings and above (suppressing info and below)
    logging.getLogger().setLevel(logging.WARNING)
    gene_anno = read_gtf(geno_anno)

    #Only include 'gene_type' == 'protein_coding' & 'feature' either 'gene' or 'exon' entries
    gene_anno = gene_anno[gene_anno['gene_type']=='protein_coding']
    gene_anno = gene_anno[gene_anno['feature'].isin(['gene', 'exon'])]

    #Remove entries where the parent "gene" contains the 'readthrough_gene' tag
    readthrough_gene_id_set = set(gene_anno[gene_anno['tag'].str.contains('readthrough_gene')]['gene_id'])
    gene_anno = gene_anno[~gene_anno['gene_id'].isin(readthrough_gene_id_set)]

    #Remove chrY and chrM entries
    gene_anno = gene_anno[(gene_anno['seqname']!='chrM') & (gene_anno['seqname']!='chrY')]

    #Make into BED format
    gene_anno['start'] = gene_anno['start'].sub(1)

    return gene_anno

def DR_analysis(geno_frag_class, dr, gene_anno_exon, output_dir, num_processes, logger):
    logger.write('-> Calculating DR scores...')
    chr_list = pd.unique(geno_frag_class['fragment_chr']).tolist()
    args_list = [(geno_frag_class[geno_frag_class['fragment_chr']==chrom], dr[dr['DR_chr']==chrom]) for chrom in chr_list]
    with multiprocessing.Pool(num_processes) as pool:
        geno_frag_class_dr = pool.starmap(get_geno_frag_class_dr, args_list)
    geno_frag_class_dr = pd.concat(geno_frag_class_dr, ignore_index=True)

    logger.write('--> Generating figure 1...')
    p_anno = []
    #eo PIF vs eo non-PIF
    p = mwu_test('eo_PIF', 'eo_non-PIF', 'fragment_class', 'fragment_DR_score', geno_frag_class_dr)
    p_anno.append(get_p_anno(p, adj=False))
    #nc PIF vs nc non-PIF
    p = mwu_test('nc_PIF', 'nc_non-PIF', 'fragment_class', 'fragment_DR_score', geno_frag_class_dr)
    p_anno.append(get_p_anno(p, adj=False))
    get_DR_score_by_fragment_class_figure1(geno_frag_class_dr, output_dir, p_anno)

    geno_frag_class_dr_closest = get_geno_frag_class_dr_closest(geno_frag_class_dr, gene_anno_exon)
    logger.write('---> Generating figure 2...')
    p_anno = []
    #nc PIF vs nc non-PIF
    p = mwu_test('nc_PIF', 'nc_non-PIF', 'fragment_class', 'distance', geno_frag_class_dr_closest)
    p_anno.append(get_p_anno(p, adj=False))
    get_distance_to_nearest_exon_figure2(geno_frag_class_dr_closest, output_dir, p_anno)

    logger.write('----> Generating figure 3...')
    #Extract distance group
    dist_group_list = pd.unique(geno_frag_class_dr_closest['distance_group']).tolist()

    #Put the Wilcoxon rank-sum test P for each distance group in the dictionary
    dist_group_anno = {}
    for dist_group in dist_group_list:
        geno_frag_class_dr_closest_dist = geno_frag_class_dr_closest[geno_frag_class_dr_closest['distance_group']== dist_group]
        dist_group_anno[dist_group] = (mwu_test('nc_PIF', 'nc_non-PIF', 'fragment_class', 'fragment_DR_score', geno_frag_class_dr_closest_dist))

    #Replace with Bonferroni corrected P for each distance group in the dictionary
    adjusted_p = multitest.multipletests(list(dist_group_anno.values()), method='bonferroni')[1]
    for i, dist_group in enumerate(dist_group_anno.keys()):
        dist_group_anno[dist_group] = adjusted_p[i]

    #Replace with P annotation for the Bonferroni corrected P
    for dist_group in dist_group_anno.keys():
        dist_group_anno[dist_group] = get_p_anno(dist_group_anno[dist_group], adj=True)
    get_DR_score_by_distance_group_figure3(geno_frag_class_dr_closest, output_dir, dist_group_anno)


#Get the DR scores for each fragment
def get_geno_frag_class_dr(geno_frag_class_chr, dr_chr):
    geno_frag_class_chr_bedtool = pybedtools.BedTool.from_dataframe(geno_frag_class_chr).sort()
    dr_chr_bedtool = pybedtools.BedTool.from_dataframe(dr_chr).sort()

    overlap_chr = dr_chr_bedtool.intersect(geno_frag_class_chr_bedtool, wa = True, wb = True, f=0.50)
    overlap_chr = overlap_chr.to_dataframe(names=dr_chr.columns.tolist() + geno_frag_class_chr.columns.tolist())
    
    if len(overlap_chr) != 0:
        #Calculate the median DR_score for each fragment...
        geno_frag_class_chr_dr = overlap_chr.groupby(['fragment_chr', 'fragment_start', 'fragment_end', 'fragment_class'])['DR_score'].median().reset_index()
        geno_frag_class_chr_dr = geno_frag_class_chr_dr.rename(columns={'DR_score':'fragment_DR_score'})

        return geno_frag_class_chr_dr
    else:
        return pd.DataFrame()

#Get the distance to the closest exon and distance group annotation for each fragment
def get_geno_frag_class_dr_closest(geno_frag_class_dr, gene_anno_exon):
    gene_anno_exon_bedtool = pybedtools.BedTool.from_dataframe(gene_anno_exon).sort()
    geno_frag_class_dr_bedtool = pybedtools.BedTool.from_dataframe(geno_frag_class_dr).sort()
    closest = geno_frag_class_dr_bedtool.closest(gene_anno_exon_bedtool, d=True)

    geno_frag_class_dr_closest = closest.to_dataframe(names=geno_frag_class_dr.columns.tolist()+gene_anno_exon.columns.tolist() + ['distance'])
    geno_frag_class_dr_closest = geno_frag_class_dr_closest[['fragment_chr', 'fragment_start', 'fragment_end', 'fragment_class', 'fragment_DR_score', 'distance']].drop_duplicates()

    bins = [x for x in range(0, 1000000, 100000)] + [1000000, float('inf')]
    geno_frag_class_dr_closest['distance_group'] = pd.cut(geno_frag_class_dr_closest['distance'], bins=bins, right=False)
    geno_frag_class_dr_closest['distance_group'] = geno_frag_class_dr_closest['distance_group'].astype('str') #Make 'distance_group' from categorical datatype into string!

    return geno_frag_class_dr_closest

#Wilcoxon rank-sum test
def mwu_test(a, b, column1, column2, df):
    a_df = df[df[column1]==a][column2]
    b_df = df[df[column1]==b][column2]
    p = mannwhitneyu(a_df, b_df, method='asymptotic')[1]
    
    return p

def get_p_anno(p, adj):
    if not adj:
        if p <= 0.0001:
            return '****'
        elif p <=0.001:
            return '***'
        elif p <= 0.01:
            return '**'
        elif p <= 0.05:
            return '*'
        else:
            return 'ns'
    else:
        if p <= 0.0001:
            return '††††'
        elif p <=0.001:
            return '†††'
        elif p <= 0.01:
            return '††'
        elif p <= 0.05:
            return '†'
        else:
            return 'ns'
        
def get_DR_score_by_fragment_class_figure1(geno_frag_class_dr, output_dir, p_anno):
    fig, ax = plt.subplots(figsize=(4,5))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.rcParams['pdf.fonttype'] = 42
    fontsize = 'large'

    order = ['promoter_fragment', 'eo_PIF', 'eo_non-PIF', 'nc_PIF', 'nc_non-PIF']
    name = ['Promoter fragment', 'eo PIF', 'eo non-PIF', 'nc PIF', 'nc non-PIF']
    colors = ['#3876b0', '#68cf51', '#d44944', '#68cf51', '#d44944']

    sns.violinplot(data=geno_frag_class_dr, x='fragment_class', y='fragment_DR_score', 
                   order=order, palette=colors, inner='box', scale='width', ax=ax, width=0.5)

    for line in ax.lines:
        line.set_color('black')

    for collection in ax.collections:
        collection.set_edgecolor('black')

    ax.tick_params(labelsize=fontsize)
    ax.set_xticklabels(labels=name, rotation=45, ha='right')

    ax.set_ylabel('Depletion rank score', fontsize=fontsize)
    ax.set_xlabel('')
    
    #Annotations for significance
    y_max = geno_frag_class_dr['fragment_DR_score'].max()
    y_significance = y_max * 1.1

    #eo PIF vs eo non-PIF
    ax.annotate("", xy=(1, y_significance), xycoords='data', xytext=(2, y_significance), textcoords='data',
                arrowprops=dict(arrowstyle="-", ec='black', connectionstyle="bar,fraction=0.2"))
    ax.text(1.5, y_significance * 1.02, p_anno[0], ha='center', va='bottom', color='black', fontsize=fontsize)

    #nc PIF vs nc non-PIF
    ax.annotate("", xy=(3, y_significance), xycoords='data', xytext=(4, y_significance), textcoords='data',
                arrowprops=dict(arrowstyle="-", ec='black', connectionstyle="bar,fraction=0.2"))
    ax.text(3.5, y_significance * 1.02, p_anno[1], ha='center', va='bottom', color='black', fontsize=fontsize)

    plt.tight_layout()

    plt.savefig(os.path.join(output_dir, 'DR_score_figure1.pdf'), dpi=300, facecolor='white')
    plt.savefig(os.path.join(output_dir, 'DR_score_figure1.png'), dpi=300, facecolor='white')

def get_distance_to_nearest_exon_figure2(geno_frag_class_dr_closest, output_dir, p_anno):
    fig, ax = plt.subplots(figsize=(2,5))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.rcParams['pdf.fonttype'] = 42
    fontsize = 'large'

    order = ['nc_PIF', 'nc_non-PIF']
    name = ['nc PIF', 'nc non-PIF']
    palette=['#68cf51', '#d44944']

    PROPS = {
        'boxprops':{'edgecolor':'black'},
        'medianprops':{'color':'black'},
        'whiskerprops':{'color':'black'},
        'capprops':{'color':'black'}}

    sns.boxplot(data=geno_frag_class_dr_closest[geno_frag_class_dr_closest['fragment_class'].isin(order)], x='fragment_class', y='distance', showfliers=False, order=order, palette=palette, ax=ax, **PROPS)

    ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(base=100000))

    ax.tick_params(labelsize=fontsize)
    ax.set_xticklabels(labels=name, rotation=45, ha='right')

    ax.set_ylabel('Distance to nearest exon (bp)', fontsize=fontsize)
    ax.set_xlabel('')
    
    # Calculate the IQR for each group
    Q1 = geno_frag_class_dr_closest.groupby('fragment_class')['distance'].quantile(0.25)
    Q3 = geno_frag_class_dr_closest.groupby('fragment_class')['distance'].quantile(0.75)
    IQR = Q3 - Q1

    # Calculate the whisker position for each group
    whisker_pos = Q3 + 1.5 * IQR

    # Find the maximum whisker position across all groups
    max_whisker_pos = whisker_pos.max()

    # Use this maximum whisker position to the annotation
    y_significance = max_whisker_pos * 1.1  # Slightly above the maximum whisker position

    # Update the ylim to make sure there's enough space for the annotation
    ax.set_ylim(ax.get_ylim()[0], y_significance * 1.1)
    
    ax.annotate("", xy=(0, y_significance), xycoords='data', xytext=(1, y_significance), textcoords='data',
                arrowprops=dict(arrowstyle="-", ec='black', connectionstyle="bar,fraction=0.2"))
    ax.text(0.5, y_significance * 1.02, p_anno[0], ha='center', va='bottom', color='black', fontsize=fontsize)

    plt.tight_layout()

    plt.savefig(os.path.join(output_dir, 'figure2.pdf'), dpi=300, facecolor='white')
    plt.savefig(os.path.join(output_dir, 'figure2.png'), dpi=300, facecolor='white')

def get_DR_score_by_distance_group_figure3(geno_frag_class_dr_closest, output_dir, dist_group_anno):
    mapping = {'[0.0, 100000.0)':'[0-100k)', 
             '[100000.0, 200000.0)':'[100k-200k)', 
             '[200000.0, 300000.0)':'[200k-300k)', 
             '[300000.0, 400000.0)':'[300k-400k)', 
             '[400000.0, 500000.0)':'[400k-500k)', 
             '[500000.0, 600000.0)':'[500k-600k)', 
             '[600000.0, 700000.0)':'[600k-700k)', 
             '[700000.0, 800000.0)':'[700k-800k)',
             '[800000.0, 900000.0)':'[800k-900k)', 
             '[900000.0, 1000000.0)':'[900k-1000k)',
             '[1000000.0, inf)':'>1000k'}

    geno_frag_class_dr_closest['distance_group'] = geno_frag_class_dr_closest['distance_group'].replace(mapping)

    fig, ax = plt.subplots(figsize=(8,6))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.rcParams['pdf.fonttype'] = 42
    fontsize = 'large'

    order = ['nc_PIF', 'nc_non-PIF']
    hue_palette=['#68cf51', '#d44944']

    sns.violinplot(data=geno_frag_class_dr_closest, x='distance_group', y='fragment_DR_score', hue='fragment_class', hue_order=order, 
                   palette=hue_palette)

    for line in ax.lines: 
            line.set_color('black')

    for collection in ax.collections:
        collection.set_edgecolor('black')

    ax.tick_params(labelsize=fontsize)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    ax.set_ylabel('Depletion rank score', fontsize=fontsize)
    ax.set_xlabel('Distance to nearest exon (bp)', fontsize=fontsize)

    a = mpatches.Patch(color=hue_palette[0], label='nc PIF')
    b = mpatches.Patch(color=hue_palette[1], label='nc non-PIF')
    plt.legend(handles=[a,b], bbox_to_anchor=(0.5, -0.5), loc='center')
    
    #Annotations for significance
    y_max = geno_frag_class_dr_closest['fragment_DR_score'].max()
    y_significance = y_max * 1.12
    
    for i, dist_group in enumerate(dist_group_anno.keys()):
        ax.annotate("", xy=(i-0.25, y_significance), xycoords='data', xytext=(i+0.25, y_significance), textcoords='data',
                    arrowprops=dict(arrowstyle="-", ec='black', connectionstyle="bar,fraction=0.2"))
        ax.text(i, y_significance * 1.04, dist_group_anno[dist_group], ha='center', va='bottom', color='black', fontsize=fontsize)

    plt.tight_layout()

    plt.savefig(os.path.join(output_dir, 'DR_score_by_distance_group_figure3.pdf'), dpi=300, facecolor='white')
    plt.savefig(os.path.join(output_dir, 'DR_score_by_distance_group_figure3.png'), dpi=300, facecolor='white')


if __name__=='__main__':
    start_time = time.time()
    pd.options.mode.chained_assignment = None

    args = parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    logger = Logger(args.output_dir)
    logger.write('SETTINGS\n========')
    for arg in vars(args):
        logger.write(f'{arg}:\t {vars(args)[arg]}')

    pybedtools.set_tempdir(os.getcwd())

    logger.write('\n1. Processing input files...')
    dr, geno_frag_class, gene_anno_exon = process_input(args.dr_ref, args.geno_frag_class, args.geno_anno)

    logger.write('\n2. Performing depletion rank analysis...')
    DR_analysis(geno_frag_class, dr, gene_anno_exon, args.output_dir, args.num_processes, logger)

    pybedtools.cleanup(remove_all=True)

    logger.write('\nDone.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60:.2f} minutes.')
    logger.close()







    



