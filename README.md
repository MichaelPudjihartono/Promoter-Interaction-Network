# Promoter Interaction Network (PIN)
  
  <img width="400" height="400" alt="Melanoma PIN circos" src="https://github.com/MichaelPudjihartono/Promoter-Interaction-Network/assets/79574643/f13b0d8f-43b3-4d63-8e95-b8dfc6c7d034">


## Characterize distal regulatory elements through mapping of genome-wide physical interactions anchored to gene promoters

## Installation
Clone the repository using the following command:

```
git clone https://github.com/MichaelPudjihartono/Promoter-Interaction-Network.git
cd Promoter-Interaction-Network
```

## Prerequisites
PIN requires Python 3.7 or higher. Ensure that your system meets this requirement before installation.

## Dependencies
PIN utilizes several Python packages for data manipulation, analysis, and visualization. Install the following packages:
- pandas: For data manipulation and analysis.
- pybedtools: A flexible Python wrapper for working with genomic data in BED format.
- gtfparse: To read and parse GTF files.
- matplotlib and seaborn: For generating plots and visualizations of data.
- scipy and statsmodels: For statistical analysis.
- pyCircos: For creating Circos plots of PIN.

## Basic Usage
```
Promoter Interation Network: Characterize distal regulatory elements through mapping of genome-wide physical interactions anchored to gene promoters.

optional arguments:
  -h, --help            show this help message and exit
  -g GENE, --gene GENE  Gene names or GENCODE IDs of anchor genes. A file with one gene per line, header should either be 'gene_name' or 'gene_id'. Default is all genes with available promoter information
                        in the promoter reference file.
  -e {MboI,DpnII,HindIII}, --enzyme {MboI,DpnII,HindIII}
                        The enzyme used to generate the fragments (e.g. HindIII)
  -hic HIC, --hic HIC   The filepath to a directory of hi-c interaction db files. Each hi-c db file should be named as 'cell-line_replicate' for the first two names (e.g. SK-
                        MEL-5_GSM2827188_merged_nodups.db), following columns are required: chr1, fragment1, chr2, fragment2.
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Directory to write results.
  -nf, --no-figures     Disables the automatic generation of circos plots and DR score analysis. By default, these analyses and plots are generated.
  -p NUM_PROCESSES, --num-processes NUM_PROCESSES
                        Number of CPUs to use (default: half the number of CPUs).
  -c CONFIG, --config CONFIG
                        The configuration file to use (default: docs/PIN.conf).
```
**How to generate Hi-C library files:** _Process your raw Hi-C data with the Juicer of HOMER pipeline (e.g. as described in PMID: 30518762) to obtain Hi-C chromatin interaction library files with the following format: read name, strand1, chr1, position1, fragment1, strand2, chr2, position2, fragment2. We recommend removing all read pairs with mapping quality score < 30. PIN works with any Hi-C pipeline that can generate interaction files in the required data format (e.g. HOMER, Juicer)._ 

## Required Reference Files
The PIN software requires several genomic reference files to function correctly. These files are too large to host on GitHub and must be generated and placed in the correct directories as described here.

Ensure you organize the reference files according to the following directory structure expected by the software:

```
lib/
├── reference_files/
    ├── promoter_reference.txt
    ├── genomic_fragment_HindIII.bed
    ├── promoter_lookup_hindiii.db
    ├── genomic_fragment_MboI.bed
    ├── promoter_lookup_mboi.db
    ├── gencode.v43.annotation.gtf.gz
    ├── DR.gor
```
| File name               | Description | How to generate |
|-------------------------|-------------|-----------------|
| `promoter_reference.txt`| This file defines the promoter regions of genes, must include the following columns: `promoter_chr`, `promoter_start`, `promoter_end`, `gene_name`, `gene_id`, `gene_strand`. Coordinates must be in BED format (0-based). | Promoter regions can be defined based on user preference. A typical approach might involve using 2,000 bp upstream and 200 bp downstream of the gene start positions as annotated in resources like GENCODE. This method aligns with the default behavior of the `promoters()` function in GenomicRanges.<br><br>In our associated publication, promoters were defined more conservatively using data from the refTSS database, which includes experimentally validated TSS sites. Users are encouraged to use this refined promoter reference for enhanced accuracy in their analyses (available as Supplementary Table 2).|
|```genomic_fragment_HindIII.bed```| A sorted file for genomic coordinates of HindIII restriction enzyme cleavage sites in BED format, must include the following columns: `fragment_chr`, `fragment_start`, `fragment_end`, `fragment`. <br><br> _Note: ```fragment``` is restriction site fragment, it must be a sorted number starting from 1 for each chromosome_.| Download the human reference genome sequence from NCBI (e.g. https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/) and in-silico digest it at HindIII restriction enzyme cleavage sites using the Restriction package of Biopython.|
|```promoter_lookup_hindiii.db```| A db file of HindIII fragments encompassing the promoter regions of genes. | A helper script (```init_promoter_fragment_db.py```) is provided to generate this file for user convenience. It needs ```promoter_reference.txt``` and ```genomic_fragment_HindIII.bed``` as inputs.|
|```genomic_fragment_MboI.bed```| Genomic coordinates of MboI restriction enzyme cleavage sites in BED format. | See how to generate ```genomic_fragment_HindIII.bed```.<br><br> _Note: DpnII has identical cleavages site to MboI, users only need to generate the MboI files if their hi-c libraries are fragmented using DpnII, the rest are handled by the software_.|
|```promoter_lookup_mboi.db```| A db file of MboI fragments encompassing the promoter regions of genes. | See how to generate ```genomic_fragment_HindIII.bed```. <br><br> _Note: DpnII has identical cleavage sites to MboI, users only need to generate the MboI files if their hi-c libraries are fragmented using DpnII, the rest are handled by the software_.|
|```gencode.v43.annotation.gtf.gz```| GENCODE genome annotation file. | Download the file from GENCODE (e.g. https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz for version 43).|
|```DR.gor```| Depletion rank score reference file. | Available as a supplementary data from Halldorsson et al. (PMID: 35859178).|

_**Note:** The file names listed above are used as examples in the provided `.conf` file. If you choose to use these example file names, no additional configuration is needed. However, if you opt for different file names (e.g. if you want a to use a different GENCODE version, or simply use shorter names), ensure you update the file paths in the `.conf` file accordingly._

## Output
PIN outputs a number of files and plots.

### Result files
| File                | Description |
|---------------------|-------------|
| `PIN.txt`           | This file is the primary output of the PIN analysis. It contains mappings of all distal genomic fragments that physically interact with gene promoters. Each line represents a distinct interaction, detailing the involved genomic fragment and its gene. |
| `PIN_gene_anno.txt` | This file augments the `PIN.txt` data by including additional genomic information about the anchor gene. Specifically, it provides the ```gene_start``` and ```gene_end``` coordinates along with the ```gene_strand``` orientation of the gene. |
| `genomic_fragment_class.txt`| This file extends the information found in the `genomic_fragment_HindIII.bed` or `genomic_fragment_MboI.bed` files by adding an additional column that describes the `fragment_class` assigned to each genomic fragment based on the PIN analysis. The `fragment_class` is a critical output that categorizes each fragment according to its interaction properties defined during the analysis. |
| `PIN.log`           | A log file generated during the software run. It records all operational messages, errors, and other diagnostic information, which can be crucial for troubleshooting and ensuring the software is functioning correctly. |



### Result Plots
| File                                    | Description |
|-----------------------------------------|-------------|
| `DR_score_figure1.png/pdf` | Violin plots displaying the distribution of depletion rank scores for each `fragment_class`. |
| `figure2.png/pdf`          | Boxplot that compares the distribution of distances from non-coding PIFs (Promoter Interacting Fragments) and coding PIFs to the nearest exon.|
| `DR_score_by_distance_group_figure3.png/pdf` | Violin plots showing the depletion rank score distribution for non-coding vs. coding PIFs across different distance groups. These plots are crucial for analyzing how distance from exons influences the depletion scores of PIFs. |
| `Circos_PIN.png`                        | A Circos plot that summarizes the distribution of `fragment class` across the genome and the interactions anchored at gene promoters. For clarity, only Trans-inter and Trans-intra (≥ 50 Mb) interactions are shown. This plot provides a comprehensive overview of genomic interactions at a macroscopic level. |

Each plot is available in PNG format for quick viewing and PDF format for high-quality output suitable for publication.




