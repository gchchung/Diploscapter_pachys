# Novel telomere detection algorithm using TideHunter
Below you will find a brief description on how to use the novel telomere detection algorithm. Example command:

```telomere_detection.py long_reads.fasta -k 4 -K 20 -l 2000 -n 40```

## Background
The motivation for this algorithm was initially to identify the telomeric repeats of _Diploscapter pachys_ and _Diploscapter coronatus_. These two species were sequenced using Nanopore and PacBio. We reasoned that if the genomic reads were not intentionally sheared, telomeres could be captured at the 5' and the 3' ends of reads. Furthermore, if _Diploscapter_ had conventional telomeres maintained by a functional telomerase, these telomeric repeats would satisfy at least three conditions usually observed for conventional telomeres:
1. Tandem repeats - they would consist of tandemly repeating sequences (due to telomerase action)
2. Stranded occupancy - clusters of telomeric repeats and their reverse complement would be found at the 3’ and 5’ ends of reads, respectively, and
3. Log-normal occupancy - the occupancy pattern of the repeats over length would resemble an inverted log-normal cumulative curve, which reflects a naturally log-normal distribution of telomere lengths.

Thus, the telomeric repeat pattern can be found by first identifying the most frequent repeat patterns at the ends of the reads, followed by the elimination of patterns that do not satisfy the 3 conditions above.

The bulk of the code was written starting 2020 and revised bit by bit through 2023. We have tested this code on long reads derived from _Diploscapter_ (Chung _& al._ 2024, in preparation), _Caenorhabditis_ ([Yoshimura _& al._ 2019, _Genome Res_](http://genome.cshlp.org/lookup/pmidlookup?view=long&pmid=31123080)) and _Meloidogyne_ ([Dai _& al._ 2023 _Nat Comm_](https://www.nature.com/articles/s41467-023-42700-w)) nematodes across both Nanopore and PacBio platforms.

## How the algorithm works
Much of our algorithm depends on the TideHunter program written by [Yan Gao, Bo Liu, Yadong Wang, and Yi Xing (Gao & al. 2019 _Bioinformatics_)](https://academic.oup.com/bioinformatics/article/35/14/i200/5529224). Our algorithm first takes a collection of long reads (from Oxford Nanopore or PacBio sequencing platforms) and partitions each read into two - the first 1000 bps and the last 1000 bps. TideHunter is run on these partitions to detect any tandem repeats of length ```k```. Next, our algorithm reads the TideHunter output and ranks the repeated sequences together with their reverse complements based on occupancy (saved as *_TideHunter_parser_OUTPUT_CONDENSED.txt). This ranking step also merges classes of repeat sequences that are circularly permuted (repeats of ATCG are considered the same class as repeats of TCGA, CGAT and GATC).

To graph the most frequently occurring repeats at the ends of long reads, our algorithm runs TideHunter again on full-length reads and parses through the output. Only the top ```n``` repeat patterns represented in the first and last 1000 bps of reads are considered. The user provides a length window ```l``` for which occupancies will be graphed - Finally, our algorithm captures the occupancy patterns in the first ```l``` bps, last ```l``` bps and the middle of all reads longer than ```2*l```, then graphs these patterns.

Finally, the user should scan the graphed occupancy patterns to determine if the stranded occupancy and log-normal occupancy can be observed for any pattern - these should be noted and tested further (eg. by FISH).

![image](https://github.com/gchchung/Diploscapter_pachys/assets/69369525/784a6498-f333-4a8a-b2f2-148f65fc4e4b)
Figure 1: scanning long genomic reads for telomeres.

## Getting started
### Prerequisites
The script is written in Python 3. Required packages include:
- [BioPython](https://biopython.org/wiki/Download) (I/O of sequence files)
- [pillow](https://pypi.org/project/Pillow/) (for graph .png output) 
- [svgutils](https://pypi.org/project/svgutils/) (for graph .svg output)
- [plotly](https://plotly.com/python/getting-started/) (for graphing)
- [kaleido](https://pypi.org/project/kaleido/) (for plotly graphs)

In addition, TideHunter v1.4.2 should be installed in the $PATH. Later versions of TideHunter may produce an ouptut whose column orders are not the same.
- [TideHunter v1.4.2](https://github.com/yangao07/TideHunter/releases)

### Usage: Graph occupancies for many repeat types
To generate occupancy graphs, covering the first and last 2000 nucleotides, for the top 40 most frequently occurring 4-mer to 20-mer terminal repeats in "long_reads.fasta"

```telomere_detection.py long_reads.fasta -k 4 -K 20 -l 2000 -n 40```

**required arguments below:**
| flag | argument |
| ---- | -------- |
| input | FASTA file containing the long reads. |
| ```-k```   | smallest repeat period (in bps) to consider |
| ```-K```   | largest repeat period (in bps) to consider |
| ```-l```   | the window for graphing the occupancy of the repeat patterns. The value of ```l``` should be determined by trial and error to fully cover the lengths of telomeres: we found that 2000 bps were sufficient for _Diploscapter_ telomeres, while _Caenorhabditis_ telomeres required at least 6000 bps. |
| ```-n```   | graph only the top ```n``` most frequently occurring repeats in the beginnings and ends of reads.|

**Output generated:**
Output files are grouped by repeat pattern length ```k``` in folders with the name ```*_k-mers```. Repeat-pattern-specific files adopt this format in their names:
```[rank of repeat]_[repeat_type]_[type of output file]_[beginning or end of reads]```

| file | example |
| ---- | ------- |
| **occupancy counts** of a specific repeat pattern and its reverse complement at the terminal ```l``` nucleotides of reads. | ```01_TTTTTT_repeat_class_head_and_tail_count.txt``` |
| **occupancy fractions** of a specific repeat pattern and its reverse complement at the terminal ```l``` nucleotides of reads. | ```01_TTTTTT_repeat_class_head_and_tail_fractions.txt``` |
| **count summary** of a specific repeat pattern and its reverse complement, at the beginning, end, and middle of reads. | ```01_TTTTTT_repeat_class_repeat_count_summary.txt``` |
| **occupancy fraction bar plot** of a specific repeat pattern and its reverse complement at the beginning, end, and middle of reads. | ```01_TTTTTT_occupancy_graph_2000nt_all_bar.png``` (or svg) |
| **occupancy fraction line plot** of a specific repeat pattern and its reverse complement at the beginning of reads. | ```01_TTTTTT_occupancy_graph_2000nt_head.png``` (or svg) |
| **occupancy fraction line plot** of a specific repeat pattern and its reverse complement at the end of reads. | ```01_TTTTTT_occupancy_graph_2000nt_head.png``` (or svg) |


To generate the occupancy graph for a single type of repeat

## References
Chung _& al._ (2024) manuscript in preparation.

Dai _& al._ (2023) _Nat Comm_. Unzipped chromosome-level genomes reveal allopolyploid nematode origin pattern as unreduced gamete hybridization.

Gao _& al._ (2019) _Bioinformatics_. TideHunter: efficient and sensitive tandem repeat detection from noisy long-reads using seed-and-chain.

Yoshimura _& al._ (2019) _Genom Res_. Recompleting the Caenorhabditis elegans genome.

