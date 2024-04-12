import pandas as pd
import pybedtools

import argparse
import os
import sys
import time
import datetime

import sqlite3

def parse_args():
    parser = argparse.ArgumentParser(description='Create an sqlite3 database of promoter fragments.')
    parser.add_argument('-f', '--fragment', required=True, help='''The filepath to the genomic fragments BED file.
                        Must have the following columns: fragment_chr, fragment_start, fragment_end, fragment.''')
    parser.add_argument('-e', '--enzyme',required=True, help='The enzyme used to generate the fragments (e.g. HindIII).')
    parser.add_argument('-p', '--promoter', required=True, help='''The filepath to a promoter reference BED file.
                        Must have the following columns: promoter_chr, promoter_start, promoter_end, gene_name, gene_id, gene_strand. 
                        If no information on gene_strand is available, fill it with '.'.''')
    parser.add_argument('-n', '--name', help='Name of db (e.g. promoter_lookup_hindiii).', default=None)
    parser.add_argument('-o', '--output-dir', help='Directory to write results.', default=os.getcwd())

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


def find_promoter_fragment(fragment, promoter, enzyme):
    fragment = pd.read_csv(fragment, sep='\t')
    fragment_bedtools = pybedtools.BedTool.from_dataframe(fragment)

    promoter = pd.read_csv(promoter, sep='\t')
    promoter = promoter[['promoter_chr', 'promoter_start', 'promoter_end', 'gene_name', 'gene_id']]
    promoter_bedtools = pybedtools.BedTool.from_dataframe(promoter)

    promoter_fragment = fragment_bedtools.intersect(promoter_bedtools, wa=True, wb=True)
    promoter_fragment = promoter_fragment.to_dataframe(names=fragment.columns.tolist() + promoter.columns.tolist())

    promoter_fragment = promoter_fragment[['fragment_chr', 'fragment', 'gene_name', 'gene_id']].drop_duplicates()
    promoter_fragment['enzyme'] = enzyme

    return promoter_fragment

def create_promoter_fragment_db(promoter_fragment, name, output_dir):
    conn = sqlite3.connect(os.path.join(output_dir, f"{name}.db")) #This create a new sqlite3 db in the specified path
    try:
        promoter_fragment.to_sql(name, conn, if_exists='replace', index=False)
        conn.commit()
    except Exception as e:
        #Handle any exception that occurs during the database saving operations
        print(f"An error occurred while saving to the database: {e}")

    finally:
        #Ensure the connection is closed even if an error occurs
        conn.close()

        #Also, save promoter_fragments to a tsv file
        promoter_fragment.to_csv(os.path.join(output_dir, f"{name}.txt"), sep='\t', index=False)

if __name__=='__main__':
    start_time = time.time()
    pd.options.mode.chained_assignment = None

    args = parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    #Check if a custom name was provided, otherwise use a default format including the enzyme name
    if args.name is None:
        args.name = f"promoter_lookup_{args.enzyme.lower()}"

    logger = Logger(args.output_dir)
    logger.write('SETTINGS\n========')
    for arg in vars(args):
        logger.write(f'{arg}:\t {vars(args)[arg]}')
    logger.write('\n')

    logger.write('Finding promoter fragments...')
    promoter_fragment = find_promoter_fragment(args.fragment, args.promoter, args.enzyme)

    logger.write('Creating promoter fragments database...')
    create_promoter_fragment_db(promoter_fragment, args.name, args.output_dir)

    logger.write('Done.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60:.2f} minutes.')
    logger.close()