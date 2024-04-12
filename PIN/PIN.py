import pandas as pd
import pybedtools

import argparse
import configparser
import os
import sys
import time
import datetime

from sqlalchemy import create_engine
from sqlalchemy.pool import NullPool
import multiprocessing

from gtfparse import read_gtf
import warnings
import logging

import depletion_rank
import circos_PIN


def parse_args():
    parser = argparse.ArgumentParser(description='''Promoter Interation Network: Characterize distal regulatory elements
                                     through mapping of genome-wide physical interactions anchored to gene promoters.''')
    
    parser.add_argument('-g', '--gene', help='''Gene names or GENCODE IDs of anchor genes. 
                        A file with one gene per line, header should either be 'gene_name' or 'gene_id'. 
                        Default is all genes with available promoter information in the promoter reference file.''', default='All')
    parser.add_argument('-e', '--enzyme',required=True, choices=['MboI', 'DpnII', 'HindIII'], help='''The enzyme used to generate the fragments (e.g. HindIII)''')
    parser.add_argument('-hic', '--hic', required=True, help='''The filepath to a directory of hi-c interaction db files.
                        Each hi-c db file should be named as 'cell-line_replicate' for the first two names (e.g. SK-MEL-5_GSM2827188_merged_nodups.db) and have the following columns: chr1, fragment1, chr2, fragment2.''')
    parser.add_argument('-o', '--output-dir', default=os.path.join(os.getcwd(), 'results'), help='''Directory to write results.''')
    parser.add_argument('-nf', '--no-figures', action='store_true', help='Disables the automatic generation of circos plots and DR score analysis. By default, these analyses and plots are generated.')
    parser.add_argument('-p', '--num-processes', type=int, default=int(multiprocessing.cpu_count()/2), help='Number of CPUs to use (default: half the number of CPUs).')
    parser.add_argument('-c', '--config', default=os.path.join(os.path.dirname(__file__), '../lib/docs/PIN.conf'), help='''The configuration file to use (default: docs/PIN.conf).''')

    return parser.parse_args()

class Logger:
    def __init__(self, output_dir):
        self.console = sys.stdout
        self.log = open(os.path.join(output_dir, 'init_promoter_fragment_db.log'), 'w')

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

class CONFIG:
    def __init__(self, config_fp):
        config = configparser.ConfigParser()
        config.read(config_fp)

        self.lib_dir = os.path.join(os.path.dirname(__file__), config.get("PATHS", "LIB_DIR"))

        self.prom_ref = os.path.join(os.path.dirname(__file__), config.get("PATHS", "PROM_REF_FP"))

        self.hindiii_geno_frag = os.path.join(os.path.dirname(__file__), config.get("PATHS", "HINDIII_GENO_FRAG_FP"))
        self.hindiii_prom_frag = os.path.join(os.path.dirname(__file__), config.get("PATHS", "HINDIII_PROM_FRAG_FP"))

        self.mboi_geno_frag = os.path.join(os.path.dirname(__file__), config.get("PATHS", "MBOI_GENO_FRAG_FP"))
        self.mboi_prom_frag = os.path.join(os.path.dirname(__file__), config.get("PATHS", "MBOI_PROM_FRAG_FP"))

        self.geno_anno = os.path.join(os.path.dirname(__file__), config.get("PATHS", "GENO_ANNO_FP"))

        self.dr_ref = os.path.join(os.path.dirname(__file__), config.get("PATHS", "DR_REF_FP"))


def parse_gene(inp_gene, prom_ref):
    #Remember that promoter reference file must have the following columns: promoter_chr, promoter_start, promoter_end, gene_name, gene_id, gene_strand.
    gene_ref = pd.read_csv(prom_ref, sep='\t')[['gene_name', 'gene_id', 'gene_strand']].drop_duplicates()

    if inp_gene == 'All':
        gene = gene_ref #Defaults to all genes
        return gene
    
    elif os.path.isfile(inp_gene):
        inp_gene = pd.read_csv(inp_gene, sep='\t')
        if 'gene_id' in inp_gene.columns:
            mode = 'gene_id'
            gene_list = pd.unique(inp_gene['gene_id']).tolist()
        elif 'gene_name' in inp_gene.columns:
            mode = 'gene_name'
            gene_list = pd.unique(inp_gene['gene_name']).tolist()
        else:
            logger.write('ERROR: Gene file must have either "gene_id" or "gene_name" as header. Exiting...')
            sys.exit(1)
        
        gene = gene_ref[gene_ref[mode].isin(gene_list)]
        missing_gene = set(gene_list) - set(gene[mode])

        # Print a warning for any missing genes
        if missing_gene:
            logger.write(f"WARNING: The following genes are not in the promoter reference file: {', '.join(missing_gene)}.\nContinuing...")
        
        return gene
    
    else:
        logger.write(f"ERROR: {inp_gene} is not a valid file path. Exiting...")
        sys.exit(1)


def find_promoter_fragment(gene, enzyme, config):
    if enzyme == 'HindIII':
        prom_frag_db = create_engine(f'sqlite:///{config.hindiii_prom_frag}', echo=False, poolclass=NullPool)
    else: #MboI and DpnII have the same promoter fragment db
        prom_frag_db = create_engine(f'sqlite:///{config.mboi_prom_frag}', echo=False, poolclass=NullPool)
    
    prom_frag_list =[]
    with prom_frag_db.connect() as con:
        for index, row in gene.iterrows():
            sql = f"SELECT gene_id, fragment_chr, fragment, enzyme FROM promoter_lookup_{enzyme.lower()} WHERE gene_id = '{row['gene_id']}'"
            prom_frag_list.append(pd.read_sql(sql, con = con))
    
    prom_frag = pd.concat(prom_frag_list)
    prom_frag = pd.merge(gene, prom_frag, how='inner', on='gene_id')

    return prom_frag

def construct_PIN(prom_frag, hic_dir, enzyme, num_processes, config, logger):
    hic_file_list = os.listdir(hic_dir)
    hic_file_len = len(hic_file_list)

    logger.write("-> Finding promoter-interacting fragments...")
    hic_file_df_group = get_hic_file_stat(hic_file_list, logger)
    logger.write(f"--> Querying {hic_file_len} hi-c files in parallel:")
    for index, row in hic_file_df_group.iterrows():
        logger.write(f"---> {row['cell_line']} ({row['cell_line_replicates']} replicates): {row['replicate']}")

    #When we find the interacting fragments of the input genes... We only want to query the connection partner of each of the UNIQUE genomic fragment classified as promoter fragment! so we drop_duplicates() as such!
    unique_prom_frag = prom_frag[['fragment_chr', 'fragment']].drop_duplicates()
    unique_prom_frag['fragment_id'] = unique_prom_frag.apply(lambda row: str(row['fragment_chr']) + ':' + str(row['fragment']), axis=1)
    prom_frag_id_list = unique_prom_frag['fragment_id'].tolist()

    #Prepare arguments for parallel processing
    args_list = [(hic_file, hic_dir, unique_prom_frag, prom_frag_id_list, enzyme) for hic_file in hic_file_list]
    with multiprocessing.Pool(num_processes) as pool:
        #When you use starmap, it waits for all the tasks to complete before returning the results. 
        interaction = pool.starmap(query_hic_file, args_list) #starmap returns a list of the results once all tasks have completed. The results are collected in a list that is returned to the caller after all function calls have been processed
    interaction = pd.concat(interaction, ignore_index=True)

    if interaction.empty:
        logger.write('WARNING: No interactions found. Exiting...')
        sys.exit(1)

    #Impute back the information on the original gene of each query unique promoter fragment
    PIN = pd.merge(interaction, prom_frag.rename(columns = {'fragment_chr':'query_chr', 'fragment':'query_fragment'}).drop(['enzyme'], axis=1), left_on = ['query_chr', 'query_fragment'], right_on = ['query_chr', 'query_fragment'], how = 'inner')
    #Filter PIN based on the number of interactions, replicates, and cell lines
    logger.write("----> Filtering promoter-interaction network...")
    PIN, geno_frag = filter_PIN(PIN, enzyme, config)

    unique_prom_frag = unique_prom_frag.drop(columns=['fragment_id']) #unique_prom_frag is also used later in determine_fragment_class()

    return PIN, unique_prom_frag, geno_frag

def get_hic_file_stat(hic_file_list, logger):
    hic_file_df = pd.DataFrame({'hic_file':hic_file_list})
    hic_file_df['cell_line'] = hic_file_df.apply(lambda row: row['hic_file'].split('_')[0], axis=1)
    hic_file_df['replicate'] = hic_file_df.apply(lambda row: row['hic_file'].split('_')[1], axis=1)

    duplicated_replicate = hic_file_df[hic_file_df['replicate'].duplicated(keep=False)]['replicate'].drop_duplicates().tolist()
    if len(duplicated_replicate) > 0:
        logger.write(f"ERROR: Duplicate hi-c replicate(s) found: {', '.join(duplicated_replicate)}. Exiting...")
        sys.exit(1)

    else:
        hic_file_df_group = hic_file_df.groupby(['cell_line']).agg({'replicate': lambda x: ', '.join(x)}).reset_index()
        hic_file_df_group['cell_line_replicates'] = hic_file_df_group.apply(lambda row: len(row['replicate'].split(', ')), axis=1)

    return hic_file_df_group

def query_hic_file(hic_file, hic_dir, unique_prom_frag, prom_frag_id_list, enzyme):
    cell_line = hic_file.split('_')[0]
    replicate = hic_file.split('_')[1]

    hic_file_db = create_engine(f"sqlite:///{os.path.join(hic_dir, hic_file)}", echo = False)
    sql = "SELECT chr1, fragment1, chr2, fragment2 FROM interactions WHERE chr1='{}' AND fragment1={} AND chr2 IN ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X')"
    rep_interaction = []

    with hic_file_db.connect() as con:
        for index, row in unique_prom_frag.iterrows():
            from_db = pd.read_sql_query(sql.format(row['fragment_chr'][3:], row['fragment']), con = con)

            if not from_db.empty:
                from_db = from_db.rename(columns={'chr1': 'query_chr', 'fragment1':'query_fragment', 'chr2': 'fragment_chr', 'fragment2': 'fragment'})
                from_db['query_chr'] = 'chr' + from_db['query_chr'].astype(str)
                from_db['fragment_chr'] = 'chr' + from_db['fragment_chr'].astype(str)
                from_db['query_fragment'] = from_db['query_fragment'].astype('int64')
                from_db['fragment'] = from_db['fragment'].astype('int64')

                rep_interaction.append(from_db)

    if len(rep_interaction) != 0:
        #Concat it
        rep_interaction = pd.concat(rep_interaction)

        #A fragment can have many interactions with the same fragment: get unique interaction pairs
        rep_interaction = rep_interaction.drop_duplicates()

        #Add new info columns
        rep_interaction['replicate'] = replicate
        rep_interaction['cell_line'] = cell_line
        rep_interaction['enzyme'] = enzyme

        #Remove promoter-promoter interactions
        #Since the query fragments are all promoter fragments, this doubles as removing self-ligating query fragments
        rep_interaction['fragment_id'] = rep_interaction.apply(lambda row: str(row['fragment_chr']) + ':' + str(row['fragment']), axis=1)
        rep_interaction = rep_interaction[~rep_interaction['fragment_id'].isin(prom_frag_id_list)]

        rep_interaction = rep_interaction.drop(columns=['fragment_id'])

        return rep_interaction
    
    else:
        return pd.DataFrame()
    
def filter_PIN(PIN, enzyme, config):
    #1. For each target_fragment-gene pair (note: gene promoter can be located in >1 fragments): count the number of interactions that support such target_fragment-gene interaction
    PIN['num_interactions'] = PIN.groupby(['fragment_chr', 'fragment', 'gene_id']).transform('size')

    #2. For each target_fragment-gene pair: count the number of replicates (of whatever cell lines) that support such target_fragment-gene interaction
    PIN['num_replicates'] = PIN.groupby(['fragment_chr', 'fragment', 'gene_id'])['replicate'].transform('nunique')

    #3. For each target_fragment-gene pair: count the number of cell lines that support such target_fragment-gene interaction
    PIN['num_cell_lines'] = PIN.groupby(['fragment_chr', 'fragment', 'gene_id'])['cell_line'].transform('nunique')

    #4. Filter PIN
    condition = ((PIN['num_interactions'] >= 2) & (PIN['num_replicates'] >= 2) & (PIN['num_cell_lines'] >=1)) #A target_fragment-promoter interaction is valid whenever it has ≥2 supporting interactions from ≥2 different replicates of ≥1 cell lines
    PIN = PIN[condition]

    #5. Impute back information on the valid target fragments
    if enzyme == 'HindIII':
        geno_frag = pd.read_csv(config.hindiii_geno_frag, sep='\t')
    else: #MboI and DpnII have the same genomic fragment
        geno_frag = pd.read_csv(config.mboi_geno_frag, sep='\t')
    PIN = pd.merge(PIN, geno_frag, how='inner', on=['fragment_chr', 'fragment'])

    #6. Include only appropriate columns and drop_duplicates()
    PIN = PIN[['gene_name', 'gene_id', 'fragment_chr', 'fragment_start', 'fragment_end', 'fragment', 'num_interactions', 'num_replicates', 'num_cell_lines']]
    PIN = PIN.drop_duplicates()

    #7. Sort PIN
    chr_order = ['chr'+str(i) for i in range(1,23)] + ['chrX']
    chr_cat = pd.CategoricalDtype(categories=chr_order, ordered=True)
    PIN['fragment_chr'] = PIN['fragment_chr'].astype(chr_cat)
    PIN = PIN.sort_values(by=['gene_name', 'fragment_chr', 'fragment_start'])

    return PIN, geno_frag #geno_frag is also used later in determine_fragment_class()

def determine_genomic_fragment_class(PIN, unique_prom_frag, geno_frag, geno_anno, logger):
    logger.write("-> Processing genome annotation file...")
    gene_anno = process_genome_annotation(geno_anno)

    #Define promoter-interacting fragments
    PIF = PIN[['fragment_chr','fragment']].drop_duplicates()

    #Define exon fragments
    gene_anno_exon = gene_anno[gene_anno['feature']=='exon'][['seqname', 'start', 'end']].drop_duplicates()
    gene_anno_exon_bedtool = pybedtools.BedTool.from_dataframe(gene_anno_exon).sort()

    geno_frag_bedtool = pybedtools.BedTool.from_dataframe(geno_frag).sort()

    exon_overlap = gene_anno_exon_bedtool.intersect(geno_frag_bedtool, wa=True, wb=True)
    exon_overlap = exon_overlap.to_dataframe(names = gene_anno_exon.columns.tolist() + geno_frag.columns.tolist())
    exon_frag = exon_overlap[['fragment_chr', 'fragment']].drop_duplicates()

    logger.write("--> Annotating genomic fragments with fragment class...")
    #Annotate geno_frag_class by annotating geno_frag with the following fragment classes:
    #1. Annotation -> promoter_fragment
    geno_frag_class = pd.merge(geno_frag, unique_prom_frag, on=['fragment_chr', 'fragment'], how='outer', indicator = True)
    geno_frag_class['promoter_fragment'] = geno_frag_class.apply(lambda row: True if row['_merge'] == 'both' else False, axis=1)
    geno_frag_class = geno_frag_class.drop(columns=['_merge'])

    #2. Annotation -> PIF
    geno_frag_class = pd.merge(geno_frag_class, PIF, on=['fragment_chr', 'fragment'], how='outer', indicator = True) 
    geno_frag_class['PIF'] = geno_frag_class.apply(lambda row: True if row['_merge'] == 'both' else False, axis=1)
    geno_frag_class = geno_frag_class.drop(columns=['_merge'])

    #3. Annotation -> exon_fragment
    geno_frag_class = pd.merge(geno_frag_class, exon_frag, on=['fragment_chr', 'fragment'], how='outer', indicator = True)
    geno_frag_class['exon_fragment'] = geno_frag_class.apply(lambda row: True if row['_merge'] == 'both' else False, axis=1)
    geno_frag_class = geno_frag_class.drop(columns=['_merge'])

    #Determine fragment class for each genomic fragment
    geno_frag_class['fragment_class'] = geno_frag_class.apply(lambda row: 'promoter_fragment' if row['promoter_fragment'] == True else
                                                              'eo_PIF' if row['promoter_fragment']== False and row['PIF'] == True and row['exon_fragment'] == True else
                                                              'nc_PIF' if row['promoter_fragment']== False and row['PIF'] == True and row['exon_fragment'] == False else
                                                              'eo_non-PIF' if row['promoter_fragment']== False and row['PIF'] == False and row['exon_fragment'] == True else
                                                              'nc_non-PIF' if row['promoter_fragment']== False and row['PIF'] == False and row['exon_fragment'] == False else
                                                              'unannotated', axis = 1)

    #Include relevant columns only
    geno_frag_class = geno_frag_class[['fragment_chr', 'fragment_start', 'fragment_end', 'fragment', 'fragment_class']]

    return geno_frag_class, gene_anno, gene_anno_exon 
    #gene_anno is also used later in annotate_PIN()
    #gene_anno_exon is also used later in depletion_rank.DR_analysis()


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

def annotate_PIN(PIN, geno_frag_class, gene_anno, logger):
    logger.write("-> Annotating PIN with genomic fragment class...")
    PIN = pd.merge(PIN, geno_frag_class[['fragment_chr', 'fragment', 'fragment_class']], how='left')

    logger.write("--> Annotating PIN with interaction distance...")
    gene_anno_gene = gene_anno[gene_anno['feature']=='gene'][['gene_id', 'strand', 'seqname', 'start', 'end']].drop_duplicates()
    PIN = pd.merge(PIN, gene_anno_gene.rename(columns={'strand':'gene_strand', 'seqname':'gene_chr', 'start':'gene_start', 'end':'gene_end'}), how='left', on='gene_id')
    PIN['distance'] = PIN.apply(lambda row: get_PIN_distance(row), axis=1)

    logger.write("---> Annotating PIN with interaction type...")
    PIN['interaction_type'] = PIN.apply(lambda row: get_PIN_interaction_type(row), axis=1)

    PIN_gene_anno = PIN[['gene_name', 'gene_chr', 'gene_start', 'gene_end', 'gene_id', 'gene_strand', 'fragment_chr', 'fragment_start', 'fragment_end', 'fragment', 'fragment_class', 'distance', 'interaction_type', 'num_interactions', 'num_replicates', 'num_cell_lines']]
    #Include only appropriate columns!
    PIN = PIN[['gene_name', 'gene_id', 'fragment_chr', 'fragment_start', 'fragment_end', 'fragment', 'fragment_class', 'distance', 'interaction_type', 'num_interactions', 'num_replicates', 'num_cell_lines']]

    return PIN, PIN_gene_anno #PIN_gene_anno is used later in generate_circos_plot()

def get_PIN_distance(row):
    #We have to handle this in a special way since we are using BED format in our dataframes, meaning that 'start' is 0-based and 'end' exclusive (in 0-based lense)
    #For this purposes, i'll measure distance as-if they are all 0-based!! 
    #Thus every 'end' coordinate must be -1 (The same distance measurement can be achieved if we consistently use 1-based coordinate meaning that every 'start' coordinate must be +1, but in this case i'll use 0-based coordinate)
    distance = ''
    if row['gene_strand'] == '+':
        gene_start_0_based_coord = row['gene_start']
    elif row['gene_strand'] == '-':
        gene_start_0_based_coord = row['gene_end'] - 1
        
    if str(row['gene_chr']) == str(row['fragment_chr']):
        if (gene_start_0_based_coord >= row['fragment_start']) and (gene_start_0_based_coord <= (row['fragment_end']-1)):
            distance = 0
        else:
            distance = min(abs(int(gene_start_0_based_coord-row['fragment_start'])), 
                      abs(int(gene_start_0_based_coord-(row['fragment_end']-1)))
                      ) 
    else:
        distance = 'NA'
        
    return distance

def get_PIN_interaction_type(row):
    try:
        if int(row['distance']) < 1000000: #The fact that we do int() here make it appropriate to use ValueError as exception since int(string) can't be assessed using numbers operator like '<' 
            return 'Cis'
        else:
            return 'Trans-intrachromosomal'
    except ValueError: # Interchromosomal interactions have distance='NA'
        return 'Trans-interchromosomal'

if __name__=='__main__':
    start_time = time.time()
    pd.options.mode.chained_assignment = None

    args = parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    config = CONFIG(args.config)
    
    logger = Logger(args.output_dir)
    logger.write('SETTINGS\n========')
    for arg in vars(args):
        logger.write(f'{arg}:\t {vars(args)[arg]}')

    pybedtools.set_tempdir(os.getcwd())

    logger.write('\n1. Parsing gene file...')
    gene = parse_gene(args.gene, config.prom_ref)

    logger.write('\n2. Finding promoter fragments...')
    prom_frag = find_promoter_fragment(gene, args.enzyme, config)

    logger.write('\n3. Constructing promoter-interaction network...')
    PIN, unique_prom_frag, geno_frag = construct_PIN(prom_frag, args.hic, args.enzyme, args.num_processes, config, logger)

    logger.write('\n4. Determining genomic fragment classes...')
    geno_frag_class, gene_anno, gene_anno_exon = determine_genomic_fragment_class(PIN, unique_prom_frag, geno_frag, config.geno_anno, logger)
    geno_frag_class.drop(columns=['fragment']).to_csv(os.path.join(args.output_dir, 'genomic_fragment_class.txt'), sep='\t', index=False)

    logger.write('\n5. Annotating promoter-interaction network...')
    PIN, PIN_gene_anno = annotate_PIN(PIN, geno_frag_class, gene_anno, logger)
    PIN.drop(columns=['fragment']).to_csv(os.path.join(args.output_dir, 'PIN.txt'), sep='\t', index=False)
    PIN_gene_anno.drop(columns=['fragment']).to_csv(os.path.join(args.output_dir, 'PIN_gene_anno.txt'), sep='\t', index=False)

    if not args.no_figures:
        logger.write('\n6. Performing depletion rank analysis...')
        dr = pd.read_csv(config.dr_ref, sep='\t')
        dr = dr.rename(columns={'Chr':'DR_chr', 'Fromx':'DR_start', 'To':'DR_end', 'rank':'DR_score'})
        dr['DR_start'] = dr['DR_start'].sub(1)
        depletion_rank.DR_analysis(geno_frag_class, dr, gene_anno_exon, args.output_dir, args.num_processes, logger)

        logger.write('\n7. Generating circos plot of PIN...')
        circos_PIN.get_circos_plot(geno_frag_class, PIN_gene_anno, args.output_dir, logger)

    pybedtools.cleanup(remove_all=True)

    logger.write('\nDone.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60:.2f} minutes.')
    logger.close()