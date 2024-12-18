# TeloSearch-LR: a novel telomere detection algorithm using long sequencing reads without assembly
Below you will find a brief description on how to use this novel telomere detection algorithm. Example command:

```python3 TeloSearchLR.py -f long_reads.fasta -k 4 -K 20 -t 1000 -m 1 -M 100 -n 6000```

## Background
This algorithm was initially created to identify the telomeric repeat motifs of _Diploscapter pachys_ and _Diploscapter coronatus_, which have divergent telomeric sequences different from related nematodes. With the availability of Nanopore and PacBio reads for these two species, we reasoned that if the genomic reads were not intentionally sheared, telomeres should be captured at the 5' and the 3' ends of some reads. Furthermore, if _Diploscapter_ had conventional telomeres maintained by a functional telomerase, these telomeric repeats would have 3 properties usually observed for conventional telomeres:
1. Tandem repeats - they would consist of tandemly repeating sequences (due to telomerase action)
2. Stranded occupancy - clusters of telomeric repeats will be found almost exclusively at the 3' ends of reads, while clusters of their reverse complement would be found almost exclusiely at the 5’ ends of reads, and
3. Log-normal occupancy - the occupancy pattern of the repeats over length would resemble an inverted log-normal cumulative curve, which reflects a naturally log-normal distribution of telomere lengths.

Thus, the telomeric repeat pattern can be found by first identifying the most frequent repeat patterns at the ends of the reads, followed by the elimination of patterns that do not satisfy the 3 conditions above.

The bulk of the code was written starting 2020 and revised bit by bit through 2023. We have tested this code on long reads derived from _Diploscapter_ (Chung _& al._ 2024, in preparation), and the following libraries deposited in the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra).

| SRA run accesssion | species | reference | sequencing platform |
| ---- | -------- | ------------ | ---------- |
| SRR7594465 | _Caenorhabditis elegans_ | Yoshimura _& al. Genome Res 29,_ 1009–1022 (2019) | PacBio |
| SRR7874309 | _Candida albicans_ | Panthee _& al. Genome Biol Evol 10,_ 2017–2022 (2018) | ONT |
| ERR10466726 | _Kluyveromyces lactis_ | | ONT |
| ERR3313355 | _Vespa velutina_ | Darwin Tree of Life Project Consortium. _Proc Natl Acad Sci 119,_ e2115642118 (2022) | PacBio |
| SRR12427674 | _Saccharomyces cerevisiae_ | Barbitoff _& al. G3 11,_ jkab029 (2021) | ONT |
| SRR23362189 | _Magnusiomyces capitatus_ | | PacBio |
| ERR6907800 | _Diabrotica virgifera_ | [Pest Genomics Initiative - Diabrotica virgifera](https://www.pestgenomics.org/species/diabrotica-virgifera) | PacBio |
| ERR10812864 | _Poecilus cupreus_ | Darwin Tree of Life Project Consortium. _Proc Natl Acad Sci 119,_ e2115642118 (2022) | PacBio |
| SRR9852240 to SRR9852244 | _Sitophilus oryzae_ | Parisot _& al. BMC Biol._ 19, 241 (2021) | ONT |
| ERR11673226 | _Geotrupes spiniger_ | Darwin Tree of Life Project Consortium. _Proc Natl Acad Sci 119,_ e2115642118 (2022) | PacBio |
| ERR11843429 | _Hyposter dolosus_ | Darwin Tree of Life Project Consortium. _Proc Natl Acad Sci 119,_ e2115642118 (2022) | PacBio |
| SRR23734611 | _Capsicum chinense_ | Liu _& al. Nat Commun 14,_ 5487 (2023). | PacBio |
| SRR26636576 | _Sorghum bicolor_ | | ONT (ultralong) |
| SRR25177361 and SRR25177362  | _Strongyloides sterocoralis_ | | ONT |


## How the algorithm works
Much of our algorithm depends on the TideHunter program written by Yan Gao, Bo Liu, Yadong Wang, and Yi Xing ([Gao & al. 2019 _Bioinformatics_](https://academic.oup.com/bioinformatics/article/35/14/i200/5529224)). Our algorithm first takes a collection of long reads (from Oxford Nanopore or PacBio sequencing platforms) and partitions each read into two - the first 1000 bps and the last 1000 bps. TideHunter is run on these partitions to detect any tandem repeats of length ```k```. Next, our algorithm reads the TideHunter output and ranks the repeated sequences together with their reverse complements based on occupancy (saved as *_TideHunter_parser_OUTPUT_CONDENSED.txt). This ranking step also merges classes of repeat sequences that are circularly permuted (repeats of ATCG are considered the same class as repeats of TCGA, CGAT and GATC).

To graph the most frequently occurring repeats at the ends of long reads, our algorithm runs TideHunter again on full-length reads and parses through the output. Only the ```m```-th most frequent to ```M```-th most frequent repeat patterns in the first and last 1000 bps of reads are considered. The user provides a length window ```n``` for which occupancies will be graphed. Finally, our algorithm graphs the occupancy patterns in the first ```n``` bps, last ```n``` bps, and the middle of all reads longer than ```2*n```.

The user should then visually scan the graphed occupancy patterns to determine if the stranded occupancy can be observed for any repeat - these should be noted and tested further (eg. by FISH) to verify that they are telomeric repeat motifs.

![image](https://github.com/gchchung/Diploscapter_pachys/assets/69369525/8366e17a-348b-4b44-b584-380ae60c4c9a)

**Figure 1:** scanning long genomic reads for telomeres.

## Getting started
### Prerequisites
#### Manual installation of required packages
The script is written in Python 3. Required packages include the following:
- [BioPython](https://biopython.org/wiki/Download) (I/O of sequence files)
- [pillow](https://pypi.org/project/Pillow/) (for graph .png output) 
- [svgutils](https://pypi.org/project/svgutils/) (for graph .svg output)
- [plotly](https://plotly.com/python/getting-started/) (for graphing)
- [kaleido](https://pypi.org/project/kaleido/) (for plotly graphs)

In addition, TideHunter v1.4.2 should be installed in the $PATH. Later versions of TideHunter may produce an ouptut whose column orders are not the same and will be incorrectly parsed by our algorithm.
- [TideHunter v1.4.2](https://github.com/yangao07/TideHunter/releases)

Finally, the read files you will use may not already be in FASTA format, but rather in FASTQ. To convert the FASTQ file into FASTA, reformat the file and discard the quality information by running ```sed``` in your shell. For example:

```sed -n '1~4s/^@/>/p;2~4p' long_reads.fastq > long_reads.fasta```

### Installation
Required Python packages can be installed with your favourite package manager (eg. ```pip3 install [package name]``` or ```conda install [package name]```). TideHunter can be installed by downloading the source files and compiling following the developers' instructions, or by downloading a pre-compiled binary if your system is x64 Linux. To run the telomere detection algorithm, download telomere_detection.py from here and add the path to $PATH, or place the Python script directly in the folder with the reads to be analysed.

### Hardware requirements
The code was written and tested on a personal computer with Intel Core i5-7300HQ with 32 GB of RAM running Ubuntu in WSL. We expect any moderately modern computer running Linux or Unix-like environments will be sufficient.

## Usage

### Scenario 1: Graph occupancies for many repeat types to scan for candidate telomeric repeats
To generate occupancy graphs covering the first and last 6000 nucleotides of all reads, for the top 40 most frequently occurring 4-mer to 20-mer terminal repeats in "long_reads.fasta"

```python3 telomere_detection.py long_reads.fasta -k 4 -K 20 -n 6000 -m 1 -M 40```

**required arguments below:**
| flag | argument |
| ---- | -------- |
| ```-f``` | FASTA file containing the long reads. |
| ```-k```   | smallest repeat period (in bps) to consider. |
| ```-K```   | largest repeat period (in bps) to consider. |
| ```-t```   | the terminal region (in bps) to consider while ranking repeat motifs. Default = 1000 |
| ```-n```   | the window (in bps) for graphing the occupancy of the repeat patterns. The value of ```n``` should be determined by trial and error to fully cover the lengths of telomeres. |
| ```-m``` and ```-M```   | graph the occupancy patterns from the m-th most frequent motif to the M-th most frequent motif. |

**Output generated:**
Intermediate files (TideHunter outputs) are saved in the current folder. Final output files are further organized by repeat pattern length ```k``` in folders with the name ```*_k-mers```. Repeat-pattern-specific files adopt this format in their names:
```[rank of repeat]_[repeat_type]_[type of output file]_[beginning or end of reads]```



| file | example |
| ---- | ------- |
| **occupancy counts** by nucleotide position of a specific repeat pattern and its reverse complement at the terminal ```n``` nucleotides of reads. | ```01_TTTTTT_repeat_class_head_and_tail_count.txt``` |
| **occupancy fractions** by nucleotide position of a specific repeat pattern and its reverse complement at the terminal ```n``` nucleotides of reads. | ```01_TTTTTT_repeat_class_head_and_tail_fractions.txt``` |
| **count summary** of a specific repeat pattern and its reverse complement, at the beginning, end, and middle of reads. Includes a count of discarded reads (reads that are not ```n``` bps long). | ```01_TTTTTT_repeat_class_repeat_count_summary.txt``` |
| **occupancy fraction bar plot** of a specific repeat pattern and its reverse complement at the beginning, end, and middle of reads. | ```01_TTTTTT_occupancy_graph_2000nt_all_bar.png``` (or svg) |
| **occupancy fraction line plot** by nucleotide position of a specific repeat pattern and its reverse complement at the beginning of reads. | ```01_TTTTTT_occupancy_graph_2000nt_head.png``` (or svg) |
| **occupancy fraction line plot** by nucleotide position of a specific repeat pattern and its reverse complement at the end of reads. | ```01_TTTTTT_occupancy_graph_2000nt_head.png``` (or svg) |
| **occupancy and occurrence** of specific repeat patterns in the first and last 1000 bps of reads. | ```*_TideHunter_parser_OUTPUT.txt``` |
| **ranked occupancy and occurence** of specific repeat patterns and their reverse complements in the first and last 1000 bps of reads. | ```*_TideHunter_parser_OUTPUT_CONDENSED.txt``` |

Finally, a collage of all the plots with the same k-value is saved as ```*_patterns_collage.png``` (or .svg) in individual k-mer folders and a collage of all plots is saved as ```all_patterns_collage.png``` (or .svg) in the current folder. The collages show the patterns by repeat occupancy rank (top-ranked at the top, r-th ranked at the bottom), and repeat period (```k``` at the left, ```K``` at the right)

**Sample output:**



## References
Chung _& al._ (2024) manuscript in preparation.

Dai _& al._ (2023) _Nat Comm_. Unzipped chromosome-level genomes reveal allopolyploid nematode origin pattern as unreduced gamete hybridization.

Gao _& al._ (2019) _Bioinformatics_. TideHunter: efficient and sensitive tandem repeat detection from noisy long-reads using seed-and-chain.

Yoshimura _& al._ (2019) _Genom Res_. Recompleting the Caenorhabditis elegans genome.

## Conference presentations
Results generated by this algorithm have been presented at multiple academic/industry conferences, including:
- The 23rd International C. elegans Conference (online) in 2021
- Nanopore Community Meeting 2022 (New York, USA) in 2022
- The 24th International C. elegans Conference (Glasgow, Scotland) in 2023
- EvoWorm 2024 (Vienna, Austria) in 2024

## Contact
George Chung, PhD

Research Scientist, New York University

Email: [gc95@nyu.edu](mailto:gc95@nyu.edu)

Twitter: @gchchung

Project: [https://github.com/gchchung/Diploscapter_pachys/](https://github.com/gchchung/Diploscapter_pachys/)


Last updated 2023-12-01
