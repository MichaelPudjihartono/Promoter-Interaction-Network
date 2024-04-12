import pandas as pd
import pybedtools

import argparse
import os
import sys
import time
import datetime

import pycircos
import matplotlib.pyplot as plt

import collections

def parse_args():
    parser = argparse.ArgumentParser(description='''Generate circos plot visualization of promoter-interaction network from PIN.py.''')
    
    parser.add_argument('-pin', '--pin_gene_anno', required=True, help='''Gene-annotated promoter-interaction network file from PIN.py (PIN_gene_anno.txt).''')
    parser.add_argument('-g', '--geno-frag-class', required=True, help='''Genomic fragment class file from PIN.py (genomic_fragment_class.txt).''')
    parser.add_argument('-o', '--output-dir', default=os.path.join(os.getcwd(), 'results'), help='''Directory to write results.''')
   
    return parser.parse_args()

class Logger:
    def __init__(self, output_dir):
        self.console = sys.stdout
        self.log = open(os.path.join(output_dir, 'circos_PIN.log'), 'w')

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

def get_circos_plot(geno_frag_class, PIN_gene_anno, output_dir, logger):
    chr_size = get_chr_size(geno_frag_class)
    data_link = get_data_link(PIN_gene_anno)
    chr_anno = get_chr_anno(geno_frag_class)

    #Initiate core class objects
    Garc    = pycircos.Garc
    Gcircle = pycircos.Gcircle

    logger.write('-> Drawing chromosome backbone...')
    circle = Gcircle(figsize=(8,8))
    for index, row in chr_size.iterrows():
        name = row['chrom']
        length = row['end']
        arc = Garc(arc_id=name, size=length, interspace=1, raxis_range=(890,985), labelposition=80, label_visible=True, facecolor="#FFFFFF00")
        circle.add_garc(arc)

    circle.set_garcs(0,360)
    for arc_id in circle.garc_dict:
        circle.tickplot(arc_id, raxis_range=(985,1000), tickinterval=20000000, ticklabels=None) 

    logger.write('--> Drawing genomic fragment class...')
    chr_anno['width'] = chr_anno['end'] - chr_anno['start']
    color_dict = {'promoter':'#3876b0', 'PIF':'#68cf51', 'non-PIF':'#d44944'}
    chr_anno['color'] = chr_anno.apply(lambda row: color_dict[row['fragment_class']], axis=1)

    arcdata_dict = collections.defaultdict(dict)
    for index, row in chr_anno.iterrows():
        name = row['chrom']
        start = row['start']
        width = row['width']
        color = row['color']
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"] = []
            arcdata_dict[name]["colors"] = []
        arcdata_dict[name]["positions"].append(start)
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["colors"].append(color)
    
    for key in arcdata_dict:
        circle.barplot(key, 
                       data=[1]*len(arcdata_dict[key]["positions"]), #This just generates a list of [1] with length equal to the length of the ['positions'] list inside the dictionary of name #e.g. [1,1,1,1,1,1,1,1,1]
                       positions=arcdata_dict[key]["positions"], #This just extract the ['positions'] list (aka the start coordinates) inside the dictionary of name
                       width=arcdata_dict[key]["widths"], #This just extract the ['widths'] list (aka the lengths) inside the dictionary of name
                       raxis_range=[890,985], 
                       facecolor=arcdata_dict[key]["colors"])#This just extract the ['colors'] list (aka the the color of cytobands) inside the dictionary of name
    
    logger.write('---> Drawing promoter interactions...')
    trans_inter_link = data_link[data_link['interaction_type']=='Trans-interchromosomal']
    trans_inter_link = trans_inter_link.iloc[:, :6]
    values_all = [] 
    arcdata_dict = collections.defaultdict(dict)
    for index, row in trans_inter_link.iterrows():
        name1 = row['gene_chr']
        start1 = row['gene_start']
        end1 = row['gene_end']
        name2 = row['fragment_chr']
        start2 = row['fragment_start']
        end2 = row['fragment_end']
        source = (name1, start1, end1, 890)
        destination = (name2, start2, end2, 890)
        circle.chord_plot(source, destination, facecolor='#7471B3')

    trans_intra_link = data_link[data_link['interaction_type']=='Trans-intrachromosomal']
    trans_intra_link = trans_intra_link[trans_intra_link['distance'] >= 50000000] #Only draw trans-intra >= 50Mb
    trans_intra_link = trans_intra_link.iloc[:, :6]
    values_all = []
    arcdata_dict = collections.defaultdict(dict)
    for index, row in trans_intra_link.iterrows():
        name1 = row['gene_chr']
        start1 = row['gene_start']
        end1 = row['gene_end']
        name2 = row['fragment_chr']
        start2 = row['fragment_start']
        end2 = row['fragment_end']
        source = (name1, start1, end1, 890)
        destination = (name2, start2, end2, 890)
        circle.chord_plot(source, destination, facecolor='#D96013')
    
    circle.figure.savefig(os.path.join(output_dir, 'Circos_PIN.png'), dpi=300, format='png', bbox_inches='tight', facecolor='white')

def get_chr_size(geno_frag_class):
    geno_frag_class_bedtool = pybedtools.BedTool.from_dataframe(geno_frag_class)
    chr_size = geno_frag_class_bedtool.merge()
    chr_size = chr_size.to_dataframe()

    #chr_size is 1-based as per pyCircos specification (see 'example_data_chromosome_general.csv' from pyCircos docs)
    chr_size['start'] = chr_size['start'].add(1)

    return chr_size

def get_data_link(PIN_gene_anno):
    #data_link is in BED format as per pyCircos specification
    data_link = PIN_gene_anno[['gene_chr', 'gene_start', 'gene_end', 'gene_strand', 'fragment_chr', 'fragment_start', 'fragment_end', 'interaction_type', 'distance']]
    data_link['gene_start2'] = data_link.apply(lambda row: row['gene_start'] if row['gene_strand'] == '+' else row['gene_end']-65001, axis=1)
    data_link['gene_end2'] = data_link.apply(lambda row: row['gene_start2']+65000, axis=1)
    data_link = data_link[['gene_chr', 'gene_start2', 'gene_end2', 'fragment_chr', 'fragment_start', 'fragment_end', 'interaction_type', 'distance']]
    data_link = data_link.rename(columns={'gene_start2':'gene_start', 'gene_end2':'gene_end'})

    return data_link

def get_chr_anno(geno_frag_class):
    promoter = geno_frag_class[geno_frag_class['fragment_class']=='promoter_fragment']
    PIF = geno_frag_class[geno_frag_class['fragment_class'].isin(['nc_PIF', 'eo_PIF'])]
    non_PIF = geno_frag_class[geno_frag_class['fragment_class'].isin(['nc_non-PIF', 'eo_non-PIF'])]

    promoter_bedtool = pybedtools.BedTool.from_dataframe(promoter)
    PIF_bedtool = pybedtools.BedTool.from_dataframe(PIF)
    non_PIF_bedtool = pybedtools.BedTool.from_dataframe(non_PIF)

    promoter_merge = promoter_bedtool.merge()
    promoter_merge = promoter_merge.to_dataframe()
    promoter_merge['fragment_class'] = 'promoter'

    PIF_merge = PIF_bedtool.merge()
    PIF_merge = PIF_merge.to_dataframe()
    PIF_merge['fragment_class'] = 'PIF'

    non_PIF_merge = non_PIF_bedtool.merge()
    non_PIF_merge = non_PIF_merge.to_dataframe()
    non_PIF_merge['fragment_class'] = 'non-PIF'

    chr_anno = pd.concat([promoter_merge, PIF_merge, non_PIF_merge])

    chr_order = ['chr' + str(x) for x in range(1,23)] + ['chrX']
    chr_cat = pd.CategoricalDtype(categories=chr_order, ordered=True)
    chr_anno['chrom'] = chr_anno['chrom'].astype(chr_cat)
    chr_anno = chr_anno.sort_values(by=['chrom', 'start'])

    return chr_anno

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

    geno_frag_class = pd.read_csv(args.geno_frag_class, sep='\t')
    PIN_gene_anno = pd.read_csv(args.pin_gene_anno, sep='\t')

    logger.write('\n1. Generating circos plot of PIN...')
    get_circos_plot(geno_frag_class, PIN_gene_anno, args.output_dir, logger)

    logger.write('\nDone.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60:.2f} minutes.')
    logger.close()